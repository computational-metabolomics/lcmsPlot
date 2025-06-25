plot_chromatogram <- function(dataset, options, single = FALSE) {
  extra_layers <- list(
    highlight_peaks(dataset, options),
    rt_lines(options),
    legend_title(options),
    faceting(options, single),
    grid_layout(options, single)
  )
  extra_layers <- extra_layers[!sapply(extra_layers, is.null)]

  p <- ggplot(
    data = dataset$data_df,
    mapping = build_aes(x = "rt", y = "intensity", options)
  ) +
    geom_line() +
    labs(x = "RT (sec)", y = "Intensity") +
    scale_fill_discrete(guide = "none") + # Removes peak highlight legend
    theme_minimal() +
    extra_layers

  return(p)
}
