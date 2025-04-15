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

run_matching_plot_variant <- function(datasets, obj) {
  plot_config <- list(
    chromatograms = plot_chromatogram,
    mass_traces = plot_mass_trace
  )

  for (variant in plot_variants) {
    if (variant$condition(datasets, obj@options)) {
      return(variant$variant(datasets, obj, plot_config))
    }
  }

  # TODO: error handling
}

plot_data <- function(datasets, obj) {
  result_plot <- run_matching_plot_variant(datasets, obj)

  result_plot <- result_plot +
    patchwork::plot_annotation(title = obj@options$labels$title) +
    patchwork::plot_layout(axes = "collect", guides = "collect") & theme(legend.position = obj@options$legend$position)

  return(result_plot)
}


