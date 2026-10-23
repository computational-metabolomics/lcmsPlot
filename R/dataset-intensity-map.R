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
