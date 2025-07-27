plot_rt_diff <- function(datasets, dataset_type, supporting_datasets, options, single = FALSE) {
  dataset <- datasets[[dataset_type]]

  extra_layers <- list(legend_title(options))
  extra_layers <- extra_layers[!sapply(extra_layers, is.null)]

  p <- ggplot(
    data = dataset,
    mapping = build_aes(x = "rt_adj", y = "diff", options = options, group = "sample_id")
  ) +
    geom_line() +
    labs(x = "Adjusted RT (sec)", y = "Adjustment difference (sec)") +
    theme_minimal() +
    extra_layers

  return(p)
}
