plot_chromatogram <- function(datasets, dataset_type, supporting_datasets, options, single = FALSE) {
  dataset <- datasets[[dataset_type]]

  saveRDS(dataset, "test.rds")

  extra_layers <- list(
    highlight_peaks(dataset, supporting_datasets$detected_peaks, options),
    highlight_spectra_scans(datasets$spectra, options),
    rt_lines(options),
    legend_title(options),
    faceting(options, single),
    grid_layout(options, single)
  )
  extra_layers <- extra_layers[!sapply(extra_layers, is.null)]

  p <- ggplot(
    data = dataset,
    mapping = build_aes(x = "rt", y = "intensity", options)
  ) +
    geom_line() +
    labs(x = "RT (sec)", y = "Intensity") +
    scale_fill_discrete(guide = "none") + # Removes peak highlight legend
    theme_minimal() +
    extra_layers

  return(p)
}
