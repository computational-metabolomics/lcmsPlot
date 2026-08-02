#' Plot a mass trace
#'
#' `plot_mass_trace()` generates a ggplot2 from a mass trace
#'
#' @param datasets A named `list` of data frames containing the primary datasets
#' to plot. For a mass trace plot the used key is `mass_traces`.
#' @param supporting_datasets A `list` of supporting data frames, such as
#' detected peaks, used for highlighting features.
#' @param options A list of plot options, controlling units, faceting,
#' highlighting, and other visual parameters.
#' @param single A `logical` value that indicates whether it should treat
#' the plot as a single dataset variant, which can affect faceting
#' and layout behavior. Default is `FALSE`.
#' @return A `ggplot` object representing the mass trace plot.
#' @keywords internal
plot_mass_trace <- function(
    datasets,
    supporting_datasets,
    options,
    single = FALSE
) {
    dataset <- datasets$mass_traces

    extra_layers <- list(
        legend_title(options),
        faceting(options, single)
    )
    extra_layers <- remove_null_elements(extra_layers)

    p <- ggplot(
        data = dataset,
        mapping = build_aes(x = sym("rt"), y = sym("mz"), options)
    ) +
        geom_point(show.legend = single) +
        labs(x = "RT (sec)", y = "m/z") +
        theme_minimal() +
        extra_layers

    if (isTRUE(options$chrom_peak_rects$show)) {
        p <- chrom_peak_rects_layer(
            p, supporting_datasets$detected_peaks, options)
    }

    return(p)
}
