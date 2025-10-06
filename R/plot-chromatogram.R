#' Plot a chromatogram.
#'
#' @param datasets A list of data frames each corresponding to a dataset to plot.
#' @param supporting_datasets The supporting datasets (e.g., detected peaks).
#' @param options The plot object's options.
#' @param single Whether it is a single dataset variant.
#' @return The plot as a ggplot2 object.
plot_chromatogram <- function(
  datasets,
  supporting_datasets,
  options,
  single = FALSE
) {
  dataset <- datasets$chromatograms

  grouping_vars <- get_grouping_variables(options)

  dataset <- dataset %>%
    group_by(across(all_of(grouping_vars))) %>%
    mutate(
      rt_plot = case_when(
        options$chromatograms$rt_unit == "minute" ~ rt / 60,
        TRUE ~ rt),
      intensity_plot = case_when(
        options$chromatograms$intensity_unit == "relative" ~ (intensity / max(intensity)) * 100,
        TRUE ~ intensity)
    ) %>%
    ungroup()

  extra_layers <- list(
    highlight_peaks(dataset, supporting_datasets$detected_peaks, options),
    highlight_apices(dataset, options, grouping_vars),
    highlight_spectra_scans(datasets$spectra, options),
    rt_lines(options),
    legend_title(options),
    faceting(options, single),
    grid_layout(options, single)
  )
  extra_layers <- extra_layers[!sapply(extra_layers, is.null)]

  x_label <- ifelse(options$chromatograms$rt_unit == "minute", "RT (minutes)", "RT (seconds)")
  y_label <- ifelse(options$chromatograms$intensity_unit == "relative", "Relative intensity", "Intensity")

  p <- ggplot(
    data = dataset,
    mapping = build_aes(x = "rt_plot", y = "intensity_plot", options = options, group = "sample_id")
  ) +
    geom_line() +
    labs(x = x_label, y = y_label) +
    scale_fill_discrete(guide = "none") + # Removes peak highlight legend
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    theme_minimal() +
    extra_layers

  return(p)
}
