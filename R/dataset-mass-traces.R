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
