plot_spectrum <- function(datasets, dataset_type, supporting_datasets, options, single = FALSE) {
  dataset <- datasets[[dataset_type]]

  extra_layers <- list(
    legend_title(options),
    faceting(options, single)
  )
  extra_layers <- extra_layers[!sapply(extra_layers, is.null)]

  # custom_rt_labeller <- function(rt_vals) {
  #   paste0("RT: ", round(rt_vals, 3), " sec.")
  # }

  custom_rt_labeller <- function(rt_vals) {
    sapply(rt_vals, function(x) {
      x_num <- as.numeric(as.character(x))
      paste0("RT: ", round(x_num, 3), " sec.")
    })
  }

  p <- ggplot(
    data = dataset,
    mapping = build_aes(x = sym("mz"), y = sym("intensity"), options)
  ) +
    geom_segment(aes(xend = mz, yend = 0), color = "black", show.legend = single) +
    facet_wrap(~ rt, ncol = 1, labeller = labeller(rt = custom_rt_labeller)) +
    labs(x = "m/z", y = "Intensity") +
    theme_minimal() +
    extra_layers

  return(p)
}
