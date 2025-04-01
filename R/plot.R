plot_variants <- list(
  list(
    condition = function(datasets, options) length(datasets) == 1,
    variant = plot_single_dataset
  ),
  list(
    condition = function(datasets, options) is.null(options$facets$facets),
    variant = plot_multiple_datasets
  ),
  list(
    condition = function(datasets, options) !is.null(options$facets$facets),
    variant = plot_multiple_faceted_datasets
  ),
  list(
    condition = function(datasets, options) !is.null(options$grid$rows) | !is.null(options$grid$cols),
    variant = plot_multiple_gridded_datasets
  )
)

run_matching_plot_variant <- function(datasets, peakly_obj) {
  plot_config <- list(
    chromatograms = plot_chromatogram,
    mass_traces = plot_mass_trace
  )

  for (variant in plot_variants) {
    if (variant$condition(datasets, peakly_obj@options)) {
      return(variant$variant(datasets, peakly_obj, plot_config))
    }
  }

  # TODO: error handling
}

plot_data <- function(datasets, peakly_obj) {
  result_plot <- run_matching_plot_variant(datasets, peakly_obj)

  result_plot <- result_plot +
    patchwork::plot_annotation(title = peakly_obj@options$labels$title) +
    patchwork::plot_layout(axes = "collect", guides = "collect") & theme(legend.position = peakly_obj@options$legend$position)

  return(result_plot)
}


