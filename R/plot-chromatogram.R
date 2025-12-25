#' Plot chromatograms from one or more datasets
#'
#' `plot_chromatogram()` generates a ggplot2 chromatogram plot from processed
#' datasets. It can handle multiple datasets, optionally highlight detected
#' peaks or apices, and apply faceting or grid layouts based on plot options.
#'
#' @param datasets A named `list` of data frames containing the primary datasets
#' to plot. Typically includes `chromatograms` and optionally `spectra`.
#' @param supporting_datasets A `list` of supporting data frames, such as
#' detected peaks, used for highlighting features.
#' @param options A list of plot options, controlling units, faceting,
#' highlighting, and other visual parameters.
#' @param single A `logical` value that indicates whether it should treat
#' the plot as a single dataset variant, which can affect faceting
#' and layout behavior. Default is `FALSE`.
#' @return A `ggplot` object representing the chromatogram plot.
#' @keywords internal
plot_chromatogram <- function(
    datasets,
    supporting_datasets,
    options,
    single = FALSE
) {
    dataset <- datasets$chromatograms

    grouping_vars <- get_grouping_variables(options)

    dataset <- dataset |>
        group_by(across(all_of(grouping_vars))) |>
        mutate(
            rt_plot = case_when(
                options$chromatograms$rt_unit == "minute" ~ rt / 60,
                TRUE ~ rt),
            intensity_plot = case_when(
                options$chromatograms$intensity_unit == "relative" ~
                    (intensity / max(intensity)) * 100,
                TRUE ~ intensity)
        ) |>
        ungroup()

    extra_layers <- list(
        highlight_peaks(dataset, supporting_datasets$detected_peaks, options),
        highlight_apices(dataset, options, grouping_vars),
        highlight_spectra_scans(datasets$spectra, options),
        rt_lines(options),
        legend_title(options),
        faceting(options, single),
        grid_layout(options, single)
    )
    extra_layers <- remove_null_elements(extra_layers)

    x_label <- ifelse(
        options$chromatograms$rt_unit == "minute",
        "RT (minutes)",
        "RT (seconds)")
    y_label <- ifelse(
        options$chromatograms$intensity_unit == "relative",
        "Relative intensity",
        "Intensity")

    p <- ggplot(
        data = dataset,
        mapping = build_aes(
            x = "rt_plot",
            y = "intensity_plot",
            options = options,
            group = "sample_id")
    ) +
        geom_line() +
        labs(x = x_label, y = y_label) +
        scale_fill_discrete(guide = "none") + # Removes peak highlight legend
        scale_x_continuous(breaks = scales::pretty_breaks()) +
        theme_minimal() +
        extra_layers

    return(p)
}
