plot_single_dataset <- function(datasets, peakly_obj, plot_config) {
  dataset_name <- names(datasets)

  plt <- plot_config[[dataset_name]](
    data = datasets[[1]],
    options = peakly_obj@options,
    single = TRUE
  )

  plt <- patchwork::wrap_plots(plt)
  return(plt)
}

plot_multiple_datasets <- function(datasets, peakly_obj, plot_config) {
  datasets_plots <- lapply(names(datasets), function (ds_name) {
    plot_config[[ds_name]](
      data = datasets[[ds_name]],
      options = peakly_obj@options
    )
  })

  plt <- patchwork::wrap_plots(datasets_plots, ncol = 1)
  return(plt)
}

plot_multiple_faceted_datasets <- function(datasets, peakly_obj, plot_config) {
  facets_options <- peakly_obj@options$facets
  facets <- facets_options$facets

  all_plots <- peakly_obj@.data@metadata %>%
    group_by(!!!syms(facets)) %>%
    group_map(~ {
      metadata_df = .y

      datasets_subset <- lapply(datasets, function (ds) {
        ds %>% filter(across(all_of(facets), ~ . == metadata_df[[cur_column()]]))
      })

      datasets_plots <- lapply(names(datasets_subset), function (ds_name) {
        plot_config[[ds_name]](
          data = datasets_subset[[ds_name]],
          options = peakly_obj@options
        )
      })

      plot_title <- paste0(unlist(lapply(facets, function (x) metadata_df[[x]])), collapse = "\n")

      datasets_plots[[1]] <- datasets_plots[[1]] +
        ggtitle(plot_title) +
        theme(plot.title = element_text(size = 9, hjust = 0.5))

      combined_plot <- patchwork::wrap_plots(datasets_plots, ncol = 1) +
        patchwork::plot_layout(axes = "collect")

      list(
        metadata = metadata_df,
        plot = combined_plot
      )
    })

  if (is.null(facets_options$ncol) & is.null(facets_options$nrow)) {
    n_panels <- length(all_plots)
    ncol <- floor(sqrt(n_panels))
    nrow <- ceiling(n_panels / ncol)
  } else {
    ncol <- facets_options$ncol
    nrow <- facets_options$nrow
  }

  plot_grid <- lapply(all_plots, function(x) x$plot)

  plt <- patchwork::wrap_plots(plot_grid, byrow = TRUE, ncol = ncol, nrow = nrow) +
    patchwork::plot_layout(axes = "collect", guides = "collect")

  return(plt)
}

plot_multiple_gridded_datasets <- function(datasets, peakly_obj, plot_config) {
  grid_options <- peakly_obj@options$facets
  facets <- facets_options$facets

  all_plots <- peakly_obj@.data@metadata %>%
    group_by(!!!syms(facets)) %>%
    group_map(~ {
      metadata_df = .y

      datasets_subset <- lapply(datasets, function (ds) {
        ds %>% filter(across(all_of(facets), ~ . == metadata_df[[cur_column()]]))
      })

      datasets_plots <- lapply(names(datasets_subset), function (ds_name) {
        plot_config[[ds_name]](
          data = datasets_subset[[ds_name]],
          options = peakly_obj@options
        )
      })

      plot_title <- paste0(unlist(lapply(facets, function (x) metadata_df[[x]])), collapse = "\n")

      datasets_plots[[1]] <- datasets_plots[[1]] +
        ggtitle(plot_title) +
        theme(plot.title = element_text(size = 9, hjust = 0.5))

      combined_plot <- patchwork::wrap_plots(datasets_plots, ncol = 1) +
        patchwork::plot_layout(axes = "collect")

      list(
        metadata = metadata_df,
        plot = combined_plot
      )
    })

  if (is.null(facets_options$ncol) & is.null(facets_options$nrow)) {
    n_panels <- length(all_plots)
    ncol <- floor(sqrt(n_panels))
    nrow <- ceiling(n_panels / ncol)
  } else {
    ncol <- facets_options$ncol
    nrow <- facets_options$nrow
  }

  plot_grid <- lapply(all_plots, function(x) x$plot)

  plt <- patchwork::wrap_plots(plot_grid, byrow = TRUE, ncol = ncol, nrow = nrow) +
    patchwork::plot_layout(axes = "collect", guides = "collect")

  return(plt)
}
