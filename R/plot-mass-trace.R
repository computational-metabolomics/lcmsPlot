plot_mass_trace <- function(datasets, dataset_type, supporting_datasets, options, single = FALSE) {
  dataset <- datasets[[dataset_type]]

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
