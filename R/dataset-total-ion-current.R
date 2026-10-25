#' Create the total ion current dataset of a data container
#'
#' One value per MS1 scan and sample: the scan's total ion current, read from
#' the scan header or, with `options$total_ion_current$initial = FALSE`, summed
#' from the scan's peaks (see `lp_total_ion_current()`).
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
        tibble(
            intensity = .total_ion_current(
                sps, !isFALSE(options$total_ion_current$initial)),
            sample_index = s)
    }, options$parallel_param)

    .set_summary(obj, "total_ion_current", bind_rows(results))
}

# The total ion current of each scan: from the scan header when `initial` is
# TRUE, as Spectra::tic() reports it, otherwise the sum of the scan's peaks.
# Scans whose header has no value are summed from their peaks either way.
.total_ion_current <- function(sps, initial) {
    if (!initial) {
        return(Spectra::tic(sps, initial = FALSE))
    }
    values <- Spectra::tic(sps, initial = TRUE)
    missing <- is.na(values)
    if (any(missing)) {
        values[missing] <- Spectra::tic(sps[missing], initial = FALSE)
    }
    values
}
