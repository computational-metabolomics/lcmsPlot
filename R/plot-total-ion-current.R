plot_total_ion_current <- function(datasets, dataset_type, supporting_datasets, options, single = FALSE) {
  dataset <- datasets[[dataset_type]]

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
