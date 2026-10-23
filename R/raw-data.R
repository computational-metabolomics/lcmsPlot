# mzR is imported because Spectra::MsBackendMzR reads the open formats with it
# while Spectra only suggests it.

#' Load the raw spectra of a container's samples
#'
#' Raw spectra are read lazily: creating a container only standardises the
#' sample metadata, and the first layer that needs the raw data calls
#' `ensure_spectra()`. Open formats (`.mzML`, `.mzXML`, `.CDF`) are read through
#' `Spectra::MsBackendMzR`, ThermoFisher `.raw` files through
#' `MsBackendRawFileReader`. Only the scan headers are read here; peaks are
#' read on demand.
#'
#' @param container An `lcmsPlotDataContainer`.
#' @param BPPARAM A `BiocParallelParam` used to read the file headers.
#' @return The `lcmsPlotDataContainer` with the spectra linked to its samples.
#' @keywords internal
#' @importFrom mzR openMSfile
ensure_spectra <- function(container, BPPARAM = NULL) {
    if (length(MsExperiment::spectra(container@experiment)) > 0) {
        return(container)
    }

    paths <- MsExperiment::sampleData(container@experiment)$sample_path
    if (is.null(paths) || anyNA(paths)) {
        stop(
            "This layer needs the raw spectra, but the data object does not ",
            "reference raw files for every sample.")
    }

    missing_files <- paths[!file.exists(paths)]
    if (length(missing_files) > 0) {
        stop(
            "Raw file(s) not found: ",
            paste(missing_files, collapse = ", "))
    }

    BPPARAM <- .bpparam(BPPARAM)
    is_raw <- is_thermo_raw(paths)
    parts <- list()

    if (any(!is_raw)) {
        parts[[length(parts) + 1]] <- Spectra::Spectra(
            paths[!is_raw],
            source = Spectra::MsBackendMzR(),
            BPPARAM = BPPARAM)
    }

    if (any(is_raw)) {
        if (!requireNamespace("MsBackendRawFileReader", quietly = TRUE)) {
            stop(
                "Package 'MsBackendRawFileReader' is required to read the ",
                "spectra of ThermoFisher .raw files. Install it with: ",
                "BiocManager::install('MsBackendRawFileReader')")
        }
        parts[[length(parts) + 1]] <- Spectra::Spectra(
            paths[is_raw],
            source = MsBackendRawFileReader::MsBackendRawFileReader(),
            BPPARAM = BPPARAM)
    }

    sps <- if (length(parts) == 1) parts[[1]] else do.call(c, parts)
    sample_index <- match(
        normalizePath(Spectra::dataOrigin(sps), mustWork = FALSE),
        normalizePath(paths, mustWork = FALSE))

    container@experiment <- .link_spectra(
        container@experiment, sps, sample_index)
    container
}

# Attach spectra to an experiment, linking each spectrum to the sample whose
# row number is given in `sample_index`.
.link_spectra <- function(experiment, sps, sample_index) {
    if (anyNA(sample_index)) {
        stop("Some spectra could not be matched to a sample.")
    }
    MsExperiment::spectra(experiment) <- sps
    MsExperiment::linkSampleData(
        experiment,
        with = "spectra",
        sampleIndex = as.integer(sample_index),
        withIndex = seq_along(sps))
}

#' Get the raw spectra of selected samples
#'
#' @param container An `lcmsPlotDataContainer` whose spectra are loaded.
#' @param sample_index An `integer` vector of sample indices to keep, or `NULL`
#' for all samples.
#' @param ms_level An `integer` MS level to keep, or `NULL` for all levels.
#' @return A `Spectra` object with a `sample_index` spectra variable.
#' @keywords internal
raw_spectra <- function(container, sample_index = NULL, ms_level = NULL) {
    experiment <- container@experiment
    sps <- MsExperiment::spectra(experiment)
    sps$sample_index <- MsExperiment::spectraSampleIndex(experiment)

    if (!is.null(sample_index)) {
        sps <- sps[sps$sample_index %in% sample_index]
    }
    if (!is.null(ms_level)) {
        sps <- Spectra::filterMsLevel(sps, ms_level)
    }
    sps
}

#' Select the retention times to use for a set of raw spectra
#'
#' xcms keeps the adjusted retention times of an `XcmsExperiment` in the
#' `rtime_adjusted` spectra variable. Replacing `rtime` with them lets every
#' extraction work on corrected retention times without any index matching.
#'
#' @param sps A `Spectra` object.
#' @param adjusted A `logical` value indicating whether to use the adjusted
#' retention times.
#' @return The `Spectra` object, with adjusted retention times if requested.
#' @keywords internal
select_rtime <- function(sps, adjusted) {
    if (adjusted) {
        if (!"rtime_adjusted" %in% Spectra::spectraVariables(sps)) {
            stop(
                "The data object has no adjusted retention times. RT ",
                "correction requires an XcmsExperiment (or XCMSnExp) on which ",
                "adjustRtime() was run.")
        }
        Spectra::rtime(sps) <- sps$rtime_adjusted
    }
    sps
}

has_adjusted_rtime <- function(container) {
    is(container@source, "XcmsExperiment") &&
        xcms::hasAdjustedRtime(container@source)
}

.bpparam <- function(BPPARAM) {
    if (is.null(BPPARAM)) BiocParallel::SerialParam() else BPPARAM
}

# Apply `fun` to each element of `x`, in parallel when a BiocParallelParam is
# given.
.apply_samples <- function(x, fun, BPPARAM = NULL) {
    if (is.null(BPPARAM)) {
        lapply(x, fun)
    } else {
        BiocParallel::bplapply(x, fun, BPPARAM = BPPARAM)
    }
}

is_thermo_raw <- function(paths) {
    !is.na(paths) & tolower(tools::file_ext(paths)) == "raw"
}

# netCDF headers carry no usable base peak intensity, and their total ion
# current is not the sum of the scan's peaks, so both are computed from peaks.
is_netcdf <- function(paths) {
    tolower(tools::file_ext(paths)) %in% c("cdf", "nc")
}

#' Convert an `XcmsRawList` into a `Spectra` object
#'
#' Each `xcmsRaw` holds the centroids of all its scans in two vectors indexed by
#' `scanindex`. They are split per scan into an in-memory `Spectra`, with the
#' base peak intensity and total ion current of every scan precomputed.
#'
#' @param xcmsraw_list An `XcmsRawList` object.
#' @return A `Spectra` object with a `sample_index` spectra variable.
#' @keywords internal
xcmsraw_list_to_spectra <- function(xcmsraw_list) {
    parts <- lapply(seq_along(xcmsraw_list@data), function(i) {
        obj <- xcmsraw_list@data[[i]]
        mzs <- obj@env$mz
        ints <- obj@env$intensity
        n_scans <- length(obj@scantime)
        n_points <- diff(c(obj@scanindex, length(mzs)))
        scan_of_point <- factor(
            rep(seq_len(n_scans), n_points), levels = seq_len(n_scans))
        path <- if (length(obj@filepath) > 0L) obj@filepath[[1L]] else ""

        spd <- S4Vectors::DataFrame(
            msLevel = 1L,
            rtime = obj@scantime,
            scanIndex = seq_len(n_scans),
            dataOrigin = path,
            totIonCurrent = as.numeric(obj@tic),
            basePeakIntensity = vapply(
                split(ints, scan_of_point),
                function(x) if (length(x)) max(x) else 0,
                numeric(1),
                USE.NAMES = FALSE),
            sample_index = i
        )
        spd$mz <- IRanges::NumericList(
            split(mzs, scan_of_point), compress = FALSE)
        spd$intensity <- IRanges::NumericList(
            split(ints, scan_of_point), compress = FALSE)
        Spectra::Spectra(spd)
    })

    do.call(c, unname(parts))
}

# Convert an m/z range into the centre m/z and ppm tolerance that
# rawrr::readChromatogram() expects.
.mz_range_to_xic <- function(mz_range) {
    mz <- (mz_range[1] + mz_range[2]) / 2
    list(mz = mz, ppm = ((mz_range[2] - mz_range[1]) / mz) * 1e6)
}

#' Extract XICs from a ThermoFisher `.raw` file with rawrr
#'
#' Uses the vendor library's own XIC extraction, which is much faster than
#' reading the spectra. `rawrr::readChromatogram()` takes one tolerance per
#' call and formats it with `"%f"`, so windows sharing that formatted tolerance
#' are extracted together.
#'
#' @param path A `character` path to the `.raw` file.
#' @param mz_ranges A two-column `matrix` with one m/z range per row.
#' @param rt_ranges A two-column `matrix` with one RT range (seconds) per row.
#' @return A `list` with one `data.frame` (columns `rtime` and `intensity`)
#' per row of `mz_ranges`. Windows with an undefined range give an empty one.
#' @keywords internal
rawrr_xic <- function(path, mz_ranges, rt_ranges) {
    n <- nrow(mz_ranges)
    result <- rep(list(.empty_peaks()), n)
    valid <- which(
        rowSums(is.na(mz_ranges)) == 0 & rowSums(is.na(rt_ranges)) == 0)

    if (length(valid) == 0) {
        return(result)
    }

    xics <- lapply(valid, function(i) .mz_range_to_xic(mz_ranges[i, ]))
    mz <- vapply(xics, `[[`, numeric(1), "mz")
    ppm <- vapply(xics, `[[`, numeric(1), "ppm")
    groups <- split(seq_along(valid), sprintf("%f", ppm))

    for (idx in groups) {
        chroms <- rawrr::readChromatogram(
            path, type = "xic", mass = mz[idx], tol = ppm[idx][1])

        for (k in seq_along(idx)) {
            i <- valid[idx[k]]
            rt_seconds <- chroms[[k]]$times * 60
            in_range <- rt_seconds >= rt_ranges[i, 1] &
                rt_seconds <= rt_ranges[i, 2]
            result[[i]] <- data.frame(
                rtime = rt_seconds[in_range],
                intensity = chroms[[k]]$intensities[in_range])
        }
    }

    result
}

#' Read the base peak or total ion chromatogram of a `.raw` file with rawrr
#'
#' @param path A `character` path to the `.raw` file.
#' @param aggregation_fun Either `"max"` (base peak) or `"sum"` (total ion
#' current). `"mean"` is not available: rawrr reports no per-scan peak count.
#' @return A `data.frame` with columns `rtime` and `intensity`, one row per
#' MS1 scan.
#' @keywords internal
rawrr_bpc_tic <- function(path, aggregation_fun) {
    if (aggregation_fun == "mean") {
        stop(
            "aggregation_fun = \"mean\" is not available for ThermoFisher ",
            ".raw files: rawrr does not report a per-scan peak count. Use ",
            "\"max\" or \"sum\" instead.")
    }

    idx <- rawrr::readIndex(path)
    chrom <- rawrr::readChromatogram(
        rawfile = path,
        type = if (aggregation_fun == "max") "bpc" else "tic")
    ms_level <- ifelse(
        idx$MSOrder == "Ms", 1L, as.integer(sub("Ms", "", idx$MSOrder)))

    keep <- ms_level == 1L
    data.frame(
        rtime = (idx$StartTime * 60)[keep],
        intensity = chrom$intensities[keep])
}

.empty_peaks <- function() {
    data.frame(rtime = numeric(), intensity = numeric())
}
