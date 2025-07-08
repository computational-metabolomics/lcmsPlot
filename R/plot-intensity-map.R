plot_intensity_map <- function(datasets, dataset_type, supporting_datasets, options, single = FALSE) {
  dataset <- datasets[[dataset_type]]

  extra_layers <- list(
    legend_title(options),
    faceting(options, single),
    grid_layout(options, single)
  )
  extra_layers <- extra_layers[!sapply(extra_layers, is.null)]

  if (options$intensity_map$density) {
    ggplot(dataset, aes(x = .data$rt, y = .data$mz)) +
      geom_density_2d_filled(contour_var = "ndensity") +
      scale_fill_viridis_d() +
      labs(x = "RT (sec)", y = "m/z", fill = "Density") +
      theme_minimal() +
      extra_layers
  } else {
    ggplot(dataset, aes(x = .data$rt, y = .data$mz, fill = log1p(.data$intensity))) +
      geom_tile() +
      scale_fill_viridis_c() +
      labs(x = "RT (sec)", y = "m/z", fill = "log(Intensity)") +
      theme_minimal() +
      extra_layers
  }
}
