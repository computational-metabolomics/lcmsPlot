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

#' Create the mass traces of the extracted chromatograms
#'
#' A mass trace shows the raw data points (RT and m/z) inside each
#' chromatogram's extraction window. They are only built once both the
#' chromatograms and the mass trace layer are requested, and only for
#' chromatograms extracted from raw data with a finite m/z and RT window.
#'
#' @param container An `lcmsPlotDataContainer`.
#' @param options A `list` of plot options.
#' @return The `lcmsPlotDataContainer` with `summaries$mass_traces` set.
#' @keywords internal
create_mass_traces <- function(container, options) {
    chrom_data <- if (length(container@chromatograms) > 0) {
        as.data.frame(Chromatograms::chromData(container@chromatograms))
    } else {
        data.frame()
    }

    bounds <- c("mzMin", "mzMax", "rtMin", "rtMax")
    has_windows <- nrow(chrom_data) > 0 &&
        all(bounds %in% colnames(chrom_data)) &&
        !is.null(MsExperiment::sampleData(container@experiment)$sample_path)

    if (!has_windows) {
        return(.set_summary(container, "mass_traces", tibble()))
    }

    windows <- chrom_data[
        stats::complete.cases(chrom_data[, bounds]) &
            is.finite(chrom_data$mzMin) & is.finite(chrom_data$mzMax),
        intersect(
            c("sample_index", "feature_id", bounds), colnames(chrom_data)),
        drop = FALSE] |>
        distinct()

    if (nrow(windows) == 0) {
        return(.set_summary(container, "mass_traces", tibble()))
    }

    container <- ensure_spectra(container, options$parallel_param)
    points <- extract_points(
        container,
        windows,
        adjusted = has_adjusted_rtime(container),
        BPPARAM = options$parallel_param)

    .set_summary(
        container,
        "mass_traces",
        points |> select(-"intensity"))
}

#' Create the intensity maps of a data container
#'
#' An intensity map is the point cloud of raw data points (RT, m/z, intensity)
#' of each sample within an m/z and RT range, binned onto a grid unless it is
#' drawn as points.
#'
#' @param obj An `lcmsPlotDataContainer`.
#' @param options A `list` representing the plot object's options.
#' @return The `lcmsPlotDataContainer` with `summaries$intensity_maps` set: a
#' `tibble` with columns `rt`, `mz`, `intensity`, and `sample_index`.
#' @keywords internal
create_intensity_map <- function(obj, options) {
    opts <- options$intensity_maps
    samples <- .select_samples(obj, opts$sample_ids)
    obj <- ensure_spectra(obj, options$parallel_param)

    windows <- tibble(
        sample_index = samples$sample_index,
        mzMin = opts$mz_range[1],
        mzMax = opts$mz_range[2],
        rtMin = opts$rt_range[1],
        rtMax = opts$rt_range[2]
    )
    points <- extract_points(obj, windows, BPPARAM = options$parallel_param)

    if (nrow(points) > 0 && !identical(opts$geom, "point")) {
        # A scatter of the raw centroids, as xcms::plotRaw() draws, is left
        # unbinned: binning would imply continuous coverage that the
        # instrument never measured.
        points <- bin_points(points, opts$bin_rt, opts$bin_mz)
    }

    .set_summary(
        obj,
        "intensity_maps",
        points |> select("rt", "mz", "intensity", "sample_index"))
}

#' Sum the intensities of the points falling in the same grid cell
#'
#' A point cloud has millions of points and almost as many cells, which is far
#' too many groups to summarise one at a time, so the sum is grouped with
#' `data.table`.
#'
#' @param points A `data.frame` with columns `sample_index`, `rt`, `mz`, and
#' `intensity`.
#' @param bin_rt,bin_mz `numeric` values giving the RT and m/z bin widths.
#' @return A `tibble` with one row per sample and non-empty cell.
#' @keywords internal
bin_points <- function(points, bin_rt, bin_mz) {
    cells <- data.table::data.table(
        sample_index = points$sample_index,
        rt = bin_coordinate(points$rt, bin_rt),
        mz = bin_coordinate(points$mz, bin_mz),
        intensity = points$intensity)
    as_tibble(cells[
        , list(intensity = sum(intensity)),
        by = c("sample_index", "rt", "mz")])
}

#' Create the total ion current dataset of a data container
#'
#' One value per MS1 scan and sample: the scan's total ion current. It is read
#' from the scan header when the format reports it reliably, and summed from
#' the scan's peaks otherwise (netCDF files).
#'
#' @param obj An `lcmsPlotDataContainer`.
#' @param options A `list` representing the plot object's options.
#' @return The `lcmsPlotDataContainer` with `summaries$total_ion_current` set:
#' a `tibble` with columns `intensity` and `sample_index`.
#' @keywords internal
create_total_ion_current <- function(obj, options) {
    samples <- .select_samples(obj, options$total_ion_current$sample_ids)
    obj <- ensure_spectra(obj, options$parallel_param)

    results <- .apply_samples(samples$sample_index, function(s) {
        sps <- raw_spectra(obj, s, ms_level = 1L)
        tibble(intensity = .scan_intensity(sps, "sum"), sample_index = s)
    }, options$parallel_param)

    .set_summary(obj, "total_ion_current", bind_rows(results))
}

#' Create the retention time adjustment dataset of a data container
#'
#' One row per scan with its raw and adjusted retention time, from the
#' `rtime` and `rtime_adjusted` spectra variables of an `XcmsExperiment`.
#'
#' @param obj An `lcmsPlotDataContainer`.
#' @param options A `list` representing the plot object's options.
#' @return The `lcmsPlotDataContainer` with `summaries$rt_diff` set: a `tibble`
#' with columns `rt_raw`, `rt_adj`, `diff`, and `sample_index`.
#' @keywords internal
create_rt_diff <- function(obj, options) {
    sps <- raw_spectra(obj)
    rt_raw <- Spectra::rtime(sps)
    rt_adj <- sps$rtime_adjusted

    .set_summary(obj, "rt_diff", tibble(
        rt_raw = rt_raw,
        rt_adj = rt_adj,
        diff = rt_adj - rt_raw,
        sample_index = sps$sample_index
    ))
}
