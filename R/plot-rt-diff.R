#' Plot an RT alignment difference plot
#'
#' `plot_rt_diff()` generates a ggplot2 from a dataset containing the
#' differences between raw and adjusted retention times.
#'
#' @param datasets A named `list` of data frames containing the primary datasets
#' to plot. For an RT difference plot the used key is `rt_diff`.
#' @param supporting_datasets A `list` of supporting data frames, such as
#' detected peaks, used for highlighting features.
#' @param options A list of plot options, controlling units, faceting,
#' highlighting, and other visual parameters.
#' @param single A `logical` value that indicates whether it should treat
#' the plot as a single dataset variant, which can affect faceting
#' and layout behavior. Default is `FALSE`.
#' @return A `ggplot` object representing the RT difference plot.
#' @keywords internal
plot_rt_diff <- function(
    datasets,
    supporting_datasets,
    options,
    single = FALSE
) {
    dataset <- datasets$rt_diff

    extra_layers <- list(legend_title(options))
    extra_layers <- remove_null_elements(extra_layers)

    mapping <- build_aes(
        x = "rt_adj",
        y = "diff",
        options = options,
        group = "sample_id"
    )

    p <- ggplot(
        data = dataset,
        mapping = mapping
    ) +
        geom_line() +
        labs(x = "Adjusted RT (sec)", y = "Adjustment difference (sec)") +
        theme_minimal() +
        extra_layers

    return(p)
}
