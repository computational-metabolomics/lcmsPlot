#' Plot spectra.
#'
#' @param datasets A list of data frames each corresponding to a dataset to plot.
#' @param supporting_datasets The supporting datasets (e.g., detected peaks).
#' @param options The plot object's options.
#' @param single Whether it is a single dataset variant.
#' @returns The plot as a ggplot2 object.
plot_spectrum <- function(
  datasets,
  supporting_datasets,
  options,
  single = FALSE
) {
  dataset <- datasets$spectra

  extra_layers <- list(
    legend_title(options),
    faceting(options, single)
  )
  extra_layers <- extra_layers[!sapply(extra_layers, is.null)]

  dataset_for_plot <- dataset %>%
    mutate(intensity = case_when(
      reference == TRUE ~ -.data$intensity,
      TRUE ~ .data$intensity
    )) %>%
    mutate(sample_id_rt = paste0(.data$sample_id, " RT: ", round(.data$rt, 3), " sec."))

  # TODO: extract as option
  top_peaks <- dataset_for_plot %>%
    group_by(sample_id_rt) %>%
    slice_max(order_by = intensity, n = 3)

  if (length(unique(dataset_for_plot$reference)) == 2) {
    p_aes <- build_aes(x = "mz", y = "intensity", options = options, color = "reference")
  } else {
    p_aes <- build_aes(x = "mz", y = "intensity", options = options)
  }

  p <- ggplot(
    data = dataset_for_plot,
    mapping = p_aes
  ) +
    geom_segment(aes(xend = .data$mz, yend = 0), show.legend = single) +
    geom_text(
      data = top_peaks,
      aes(label = round(mz, 4)),
      nudge_y = 0.05 * max(dataset_for_plot$intensity),
      size = 3,
      color = "red"
    ) +
    facet_wrap(~ sample_id_rt, ncol = 1) +
    labs(x = "m/z", y = "Relative intensity (%)") +
    expand_limits(y = max(dataset_for_plot$intensity) * 1.1) +
    scale_y_continuous(breaks = seq(0, 100, 10)) + # TODO: extract as option
    scale_x_continuous(breaks = scales::pretty_breaks(n = 20)) + # TODO: extract as option
    theme_minimal() +
    extra_layers

  return(p)
}
