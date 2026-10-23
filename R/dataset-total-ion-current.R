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
