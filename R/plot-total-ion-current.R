#' Plot total ion current.
#'
#' @param datasets A list of data frames each corresponding to a dataset to plot.
#' @param supporting_datasets The supporting datasets (e.g., detected peaks).
#' @param options The plot object's options.
#' @param single Whether it is a single dataset variant.
#' @return The plot as a ggplot2 object.
plot_total_ion_current <- function(
  datasets,
  supporting_datasets,
  options,
  single = FALSE
) {
  dataset <- datasets$total_ion_current

  extra_layers <- list(legend_title(options))
  extra_layers <- extra_layers[!sapply(extra_layers, is.null)]

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
