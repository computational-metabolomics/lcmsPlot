is_facet <- function (o) !is.null(o$facets$facets)

is_grid <- function (o) !is.null(o$grid$rows) | !is.null(o$grid$cols)

plot_variants <- list(
  list(
    condition = function(datasets, options) length(datasets) == 1,
    variant = plot_single_dataset
  ),
  list(
    condition = function(datasets, options) !is_facet(options) & !is_grid(options),
    variant = plot_multiple_datasets
  ),
  list(
    condition = function(datasets, options) is_facet(options),
    variant = plot_multiple_faceted_datasets
  ),
  list(
    condition = function(datasets, options) is_grid(options),
    variant = plot_multiple_gridded_datasets
  )
)

run_matching_plot_variant <- function(datasets, obj) {
  plot_config <- list(
    chromatograms = plot_chromatogram,
    mass_traces = plot_mass_trace,
    spectra = plot_spectrum
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


