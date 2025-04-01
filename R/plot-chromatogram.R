plot_chromatogram <- function(data, options, single = FALSE) {
  extra_layers <- list(
    highlight_peaks(data, options),
    legend_title(options),
    faceting(options, single),
    grid_layout(options, single)
  )
  extra_layers <- extra_layers[!sapply(extra_layers, is.null)]

  p <- ggplot(
    data = data,
    mapping = build_aes(x = "rt", y = "intensity", options)
  ) +
    geom_line() +
    labs(x = "RT (sec)", y = "Intensity") +
    theme_minimal() +
    extra_layers

  return(p)
}
