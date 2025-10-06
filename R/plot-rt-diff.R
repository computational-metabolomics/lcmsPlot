#' Plot an RT alignment difference plot.
#'
#' @param datasets A list of data frames each corresponding to a dataset to plot.
#' @param supporting_datasets The supporting datasets (e.g., detected peaks).
#' @param options The plot object's options.
#' @param single Whether it is a single dataset variant.
#' @return The plot as a ggplot2 object.
plot_rt_diff <- function(
  datasets,
  supporting_datasets,
  options,
  single = FALSE
) {
  dataset <- datasets$rt_diff

  extra_layers <- list(legend_title(options))
  extra_layers <- extra_layers[!sapply(extra_layers, is.null)]

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
