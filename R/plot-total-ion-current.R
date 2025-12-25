#' Plot total ion current for one or more samples
#'
#' `plot_total_ion_current()` generates a ggplot2 of the total ion current
#' of the samples within the dataset.
#'
#' @param datasets A named `list` of data frames containing the primary datasets
#' to plot. For a TIC plot the used key is `total_ion_current`.
#' @param supporting_datasets A `list` of supporting data frames, such as
#' detected peaks, used for highlighting features.
#' @param options A list of plot options, controlling units, faceting,
#' highlighting, and other visual parameters.
#' @param single A `logical` value that indicates whether it should treat
#' the plot as a single dataset variant, which can affect faceting
#' and layout behavior. Default is `FALSE`.
#' @return A `ggplot` object representing the RT difference plot.
#' @keywords internal
plot_total_ion_current <- function(
    datasets,
    supporting_datasets,
    options,
    single = FALSE
) {
    dataset <- datasets$total_ion_current

    extra_layers <- list(legend_title(options))
    extra_layers <- remove_null_elements(extra_layers)

    geom_func_name <- paste0("geom_", options$total_ion_current$type)
    geom_func <- get(geom_func_name, asNamespace("ggplot2"))

    p <- ggplot(
        data = dataset,
        mapping = build_aes(x = "sample_id", y = "intensity", options = options)
    ) +
        geom_func() +
        theme_minimal() +
        extra_layers

    return(p)
}
