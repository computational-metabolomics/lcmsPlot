plot_spectrum <- function(datasets, dataset_type, supporting_datasets, options, single = FALSE) {
  dataset <- datasets[[dataset_type]]

  extra_layers <- list(
    legend_title(options),
    faceting(options, single)
  )
  extra_layers <- extra_layers[!sapply(extra_layers, is.null)]

  dataset_for_plot <- dataset %>%
    mutate(intensity = case_when(
      reference == TRUE ~ -intensity,
      TRUE ~ intensity
    )) %>%
    mutate(sample_id_rt = paste0(sample_id, " RT: ", round(rt, 3), " sec."))

  if (length(unique(dataset_for_plot$reference)) == 2) {
    p_aes <- build_aes(x = "mz", y = "intensity", options = options, color = "reference")
  } else {
    p_aes <- build_aes(x = "mz", y = "intensity", options = options)
  }

  p <- ggplot(
    data = dataset_for_plot,
    mapping = p_aes
  ) +
    geom_segment(aes(xend = mz, yend = 0), show.legend = single) +
    facet_wrap(~ sample_id_rt, ncol = 1) +
    labs(x = "m/z", y = "Relative intensity (%)") +
    theme_minimal() +
    extra_layers

  return(p)
}
