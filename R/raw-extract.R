#' Extract one chromatogram per window from the raw data
#'
#' The extraction engine shared by every window-based chromatogram method.
#' Each sample's raw file is read once for all of its windows: open formats
#' through `Chromatograms::ChromBackendSpectra`, which reads the peaks of the
#' scans covering the union of the windows' RT ranges in one pass, and
#' ThermoFisher `.raw` files through the vendor library's XIC extraction.
#' Samples are processed in parallel when `BPPARAM` is given. The
#' chromatograms keep the order of the windows.
#'
#' @param container An `lcmsPlotDataContainer`.
#' @param windows A `data.frame` with one row per chromatogram and columns
#' `sample_index`, `mzMin`, `mzMax`, `rtMin`, and `rtMax`, plus annotation
#' columns (such as `feature_id`) that are kept on the chromatograms.
#' @param adjusted A `logical` value indicating whether to extract on the
#' adjusted retention times.
#' @param fill_gaps A `logical` indicating whether scans with no signal in a
#' window are kept with zero intensity (`TRUE`) or dropped (`FALSE`).
#' @param ms_level An `integer` giving the MS level of the scans to use.
#' @param BPPARAM A `BiocParallelParam`, or `NULL` to process sequentially.
#' @return The `lcmsPlotDataContainer` with its `chromatograms` slot set.
#' @keywords internal
extract_chromatograms <- function(
    container,
    windows,
    adjusted = FALSE,
    fill_gaps = FALSE,
    ms_level = 1L,
    BPPARAM = NULL
) {
    windows <- as_tibble(windows)
    paths <- MsExperiment::sampleData(container@experiment)$sample_path
    use_rawrr <- .use_rawrr(paths)

    if (!all(use_rawrr[unique(windows$sample_index)])) {
        container <- ensure_spectra(container, BPPARAM)
    }

    by_sample <- split(
        seq_len(nrow(windows)),
        factor(windows$sample_index, levels = unique(windows$sample_index)))

    peaks <- .apply_samples(names(by_sample), function(sample_key) {
        rows <- by_sample[[sample_key]]
        w <- windows[rows, , drop = FALSE]
        s <- as.integer(sample_key)

        if (use_rawrr[s]) {
            rawrr_xic(
                paths[s],
                cbind(w$mzMin, w$mzMax),
                cbind(w$rtMin, w$rtMax))
        } else {
            sps <- select_rtime(
                raw_spectra(container, s, ms_level), adjusted)
            .window_peaks(sps, w, ms_level)
        }
    }, BPPARAM)

    # Back from sample order to the order of the windows.
    peaks <- unlist(peaks, recursive = FALSE, use.names = FALSE)
    peaks[unlist(by_sample, use.names = FALSE)] <- peaks
    peaks <- lapply(peaks, .apply_gap_policy, fill_gaps = fill_gaps)

    container@chromatograms <- as_memory_chromatograms(windows, peaks)
    container
}

# Use the rawrr XIC extraction for ThermoFisher files whenever rawrr is
# installed; it never needs the spectra.
.use_rawrr <- function(paths) {
    if (is.null(paths)) {
        return(logical())
    }
    is_thermo_raw(paths) & requireNamespace("rawrr", quietly = TRUE)
}

# The summed intensity per scan of every window of one sample, read in one pass
# over the sample's spectra.
.window_peaks <- function(sps, windows, ms_level) {
    if (length(sps) == 0) {
        return(rep(list(.empty_peaks()), nrow(windows)))
    }

    valid <- stats::complete.cases(
        windows[, c("mzMin", "mzMax", "rtMin", "rtMax")])
    result <- rep(list(.empty_peaks()), nrow(windows))
    if (!any(valid)) {
        return(result)
    }

    chrom_data <- data.frame(
        msLevel = as.integer(ms_level),
        dataOrigin = Spectra::dataOrigin(sps)[1],
        mzMin = windows$mzMin[valid],
        mzMax = windows$mzMax[valid],
        rtMin = windows$rtMin[valid],
        rtMax = windows$rtMax[valid])

    # Every spectrum belongs to the one sample, so it is a single chromatogram
    # group whatever the file's dataOrigin.
    Spectra::dataOrigin(sps) <- rep(chrom_data$dataOrigin[1], length(sps))
    chr <- Chromatograms::Chromatograms(
        sps, chromData = chrom_data, summarize.method = "sum")
    result[valid] <- Chromatograms::peaksData(chr)
    result
}

# Scans with no signal in the window come back with a missing intensity: keep
# them as zeros when gaps are filled, drop them otherwise.
.apply_gap_policy <- function(peaks, fill_gaps) {
    missing <- is.na(peaks$intensity)
    if (!any(missing)) {
        return(peaks)
    }
    if (fill_gaps) {
        peaks$intensity[missing] <- 0
        peaks
    } else {
        peaks[!missing, , drop = FALSE]
    }
}

#' Build an in-memory `Chromatograms` object
#'
#' @param chrom_data A `data.frame` with one row per chromatogram. Columns that
#' are entirely missing are dropped by the backend.
#' @param peaks A `list` of `data.frame`s with columns `rtime` and `intensity`,
#' one per row of `chrom_data`.
#' @return A `Chromatograms` object with a `ChromBackendMemory` backend.
#' @keywords internal
as_memory_chromatograms <- function(chrom_data, peaks) {
    if (nrow(chrom_data) == 0) {
        return(Chromatograms::Chromatograms())
    }

    chrom_data <- as.data.frame(chrom_data)
    chrom_data$sample_index <- as.integer(chrom_data$sample_index)
    if (!"msLevel" %in% colnames(chrom_data)) {
        chrom_data$msLevel <- 1L
    }
    peaks <- lapply(peaks, function(p) {
        data.frame(
            rtime = as.numeric(p$rtime),
            intensity = as.numeric(p$intensity))
    })

    Chromatograms::Chromatograms(Chromatograms::backendInitialize(
        Chromatograms::ChromBackendMemory(),
        chromData = chrom_data,
        peaksData = peaks))
}

# The base peak intensity, total ion current, or mean intensity of each scan.
.scan_intensity <- function(sps, aggregation_fun) {
    from_header <- !any(is_netcdf(Spectra::dataOrigin(sps)))

    total <- .header_value(sps, "totIonCurrent", from_header)
    if (is.null(total) && aggregation_fun != "max") {
        total <- Spectra::ionCount(sps)
        total[is.na(total)] <- 0
    }

    switch(
        aggregation_fun,
        max = {
            bpi <- .header_value(sps, "basePeakIntensity", from_header)
            if (is.null(bpi)) {
                bpi <- vapply(
                    Spectra::intensity(sps),
                    function(x) if (length(x)) max(x) else 0,
                    numeric(1))
            }
            bpi
        },
        sum = total,
        mean = {
            count <- lengths(sps)
            ifelse(count > 0, total / count, 0)
        }
    )
}

# A header column, if the format reports it for every scan. Readers mark
# unknown values as negative (mzR reports -1 for netCDF base peaks).
.header_value <- function(sps, variable, from_header) {
    if (!from_header ||
        !variable %in% Spectra::spectraVariables(sps)) {
        return(NULL)
    }
    values <- Spectra::spectraData(sps, variable)[[1]]
    if (anyNA(values) || any(values < 0)) NULL else values
}

#' Extract the raw data points inside m/z and RT windows
#'
#' Reads each sample's peaks once, for the scans covering the union of its
#' windows' RT ranges, and slices every window out of them. Shared by the mass
#' traces and the intensity maps.
#'
#' @param container An `lcmsPlotDataContainer` whose spectra are loaded.
#' @param windows A `data.frame` with one row per window and columns
#' `sample_index`, `mzMin`, `mzMax`, `rtMin`, and `rtMax`, plus annotation
#' columns (such as `feature_id`) copied onto the window's points.
#' @param adjusted A `logical` value indicating whether to use the adjusted
#' retention times.
#' @param ms_level An `integer` giving the MS level of the scans to use.
#' @param BPPARAM A `BiocParallelParam`, or `NULL` to process sequentially.
#' @return A `tibble` with columns `rt`, `mz`, `intensity`, and the window
#' columns other than the window bounds, one row per data point.
#' @keywords internal
extract_points <- function(
    container,
    windows,
    adjusted = FALSE,
    ms_level = 1L,
    BPPARAM = NULL
) {
    windows <- as_tibble(windows)
    if (nrow(windows) == 0) {
        return(tibble(rt = numeric(), mz = numeric(), intensity = numeric()))
    }

    annotation_cols <- setdiff(
        colnames(windows), c("mzMin", "mzMax", "rtMin", "rtMax"))
    by_sample <- split(
        seq_len(nrow(windows)),
        factor(windows$sample_index, levels = unique(windows$sample_index)))

    results <- .apply_samples(names(by_sample), function(sample_key) {
        w <- windows[by_sample[[sample_key]], , drop = FALSE]
        sps <- select_rtime(
            raw_spectra(container, as.integer(sample_key), ms_level),
            adjusted)
        rt <- Spectra::rtime(sps)

        in_any <- Reduce(`|`, lapply(seq_len(nrow(w)), function(k) {
            rt >= w$rtMin[k] & rt <= w$rtMax[k]
        }), FALSE)
        scans <- which(in_any)
        if (length(scans) == 0) {
            return(NULL)
        }

        scan_rt <- rt[scans]
        scan_peaks <- Spectra::peaksData(sps[scans], return.type = "list")

        bind_rows(lapply(seq_len(nrow(w)), function(k) {
            sel <- which(scan_rt >= w$rtMin[k] & scan_rt <= w$rtMax[k])
            if (length(sel) == 0) {
                return(NULL)
            }
            points <- do.call(rbind, scan_peaks[sel])
            n_points <- vapply(scan_peaks[sel], nrow, integer(1))
            keep <- points[, 1] >= w$mzMin[k] & points[, 1] <= w$mzMax[k]

            bind_cols(
                tibble(
                    rt = rep(scan_rt[sel], n_points)[keep],
                    mz = points[keep, 1],
                    intensity = points[keep, 2]),
                w[rep(k, sum(keep)), annotation_cols, drop = FALSE])
        }))
    }, BPPARAM)

    out <- bind_rows(results)
    if (nrow(out) == 0) {
        return(tibble(rt = numeric(), mz = numeric(), intensity = numeric()))
    }
    out
}

#' Read the spectra of the target scans
#'
#' @param obj An `lcmsPlotDataContainer` whose spectra are loaded.
#' @param targets A `data.frame` as returned by [spectra_targets()].
#' @param ms_level The MS level of the scans matched by retention time.
#' @param adjusted A `logical` value indicating whether to use the adjusted
#' retention times.
#' @return An in-memory `Spectra` object, one spectrum per target, with
#' spectra variables `sample_index`, `feature_id`, and `reference`. Targets
#' ordered by sample come back in their order.
#' @keywords internal
select_target_spectra <- function(obj, targets, ms_level, adjusted = FALSE) {
    if (nrow(targets) == 0) {
        return(Spectra::Spectra())
    }

    by_sample <- split(
        seq_len(nrow(targets)),
        factor(targets$sample_index, levels = unique(targets$sample_index)))

    parts <- lapply(names(by_sample), function(sample_key) {
        rows <- by_sample[[sample_key]]
        target <- targets[rows, , drop = FALSE]

        if ("scan_index" %in% colnames(target)) {
            sps <- select_rtime(
                raw_spectra(obj, as.integer(sample_key)), adjusted)
            idx <- target$scan_index
        } else {
            sps <- select_rtime(
                raw_spectra(obj, as.integer(sample_key), ms_level), adjusted)
            if (length(sps) == 0) {
                stop("No scans found with the specified MS level.")
            }
            idx <- .closest_scans(Spectra::rtime(sps), target$target_rt)
        }

        selected <- sps[idx]
        selected$feature_id <- target$feature_id
        selected
    })

    sps <- Spectra::setBackend(
        do.call(c, unname(parts)), Spectra::MsBackendMemory())
    sps$reference <- FALSE
    sps
}

# The index of the scan closest to each target retention time.
.closest_scans <- function(scan_rt, target_rt) {
    scan_order <- order(scan_rt)
    target_order <- order(target_rt)
    idx <- integer(length(target_rt))
    idx[target_order] <- scan_order[MsCoreUtils::closest(
        target_rt[target_order],
        scan_rt[scan_order],
        tolerance = Inf,
        duplicates = "keep")]
    idx
}
