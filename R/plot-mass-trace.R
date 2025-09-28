#' Plot a mass trace.
#'
#' @param datasets A list of data frames each corresponding to a dataset to plot.
#' @param supporting_datasets The supporting datasets (e.g., detected peaks).
#' @param options The plot object's options.
#' @param single Whether it is a single dataset variant.
#' @returns The plot as a ggplot2 object.
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
  extra_layers <- extra_layers[!sapply(extra_layers, is.null)]

  p <- ggplot(
    data = dataset,
    mapping = build_aes(x = sym("rt"), y = sym("mz"), options)
  ) +
    geom_point(show.legend = single) +
    labs(x = "RT (sec)", y = "m/z") +
    theme_minimal() +
    extra_layers

  return(p)
}
