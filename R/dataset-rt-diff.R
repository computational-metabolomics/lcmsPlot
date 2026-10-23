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
