plot_chromatogram <- function(datasets, dataset_type, supporting_datasets, options, single = FALSE) {
  dataset <- datasets[[dataset_type]]

  extra_layers <- list(
    highlight_peaks(dataset, supporting_datasets$detected_peaks, options),
    highlight_spectra_scans(datasets$spectra, options),
    rt_lines(options),
    legend_title(options),
    faceting(options, single),
    grid_layout(options, single)
  )
  extra_layers <- extra_layers[!sapply(extra_layers, is.null)]

  grouping_vars <- get_grouping_variables(options)
  
  dataset = dataset %>%
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
  
  highlight_apices_opts = options$chromatograms$highlight_apices
  
  top_peaks = dataset %>%
    group_by(across(all_of(grouping_vars))) %>%
    group_modify(~ {
      if (!is.null(highlight_apices_opts$column)) {
        hp <- unique(.x[[highlight_apices_opts$column]])
      } else {
        hp <- NA
      }
      
      if (!is.na(hp)) {
        # filter within ±5 rt and select max
        .x %>%
          filter(rt >= hp - 5, rt <= hp + 5) %>%
          slice_max(intensity, n = 1)
      } else if (!is.null(highlight_apices_opts$top_n)) {
        # take global maximum
        .x %>% slice_max(intensity, n = highlight_apices_opts$top_n)
      } else {
        tibble()
      }
    })
  
  x_label = ifelse(options$chromatograms$rt_unit == "minute", "RT (minutes)", "RT (seconds)")
  y_label = ifelse(options$chromatograms$intensity_unit == "relative", "Relative intensity", "Intensity")
  
  p <- ggplot(
    data = dataset,
    mapping = build_aes(x = "rt_plot", y = "intensity_plot", options = options, group = "sample_id")
  ) +
    geom_line() +
    geom_text(
      data = top_peaks,
      aes(label = round(rt_plot, 2)),
      nudge_y = 0.05 * max(dataset$intensity_plot),
      size = 3,
      color = "red"
    ) +
    labs(x = x_label, y = y_label) +
    scale_fill_discrete(guide = "none") + # Removes peak highlight legend
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    theme_minimal() +
    extra_layers
  
  return(p)
}
