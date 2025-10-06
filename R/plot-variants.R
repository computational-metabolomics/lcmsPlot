#' Plot a single dataset type.
#'
#' @param datasets A list of data frames each corresponding to a dataset to plot.
#' @param obj The lcmsPlot object.
#' @param plot_config The plot configuration.
#' A list that specifies the plot function to use for each dataset type.
#' @return The patchwork plot object.
plot_single_dataset <- function(datasets, obj, plot_config) {
  dataset_name <- names(datasets)

  plt <- plot_config[[dataset_name]](
    datasets = datasets,
    supporting_datasets = list(detected_peaks = obj@data@detected_peaks),
    options = obj@options,
    single = TRUE
  )

  plt <- patchwork::wrap_plots(plt)
  return(plt)
}

#' Plot multiple dataset types.
#'
#' @param datasets A list of data frames each corresponding to a dataset to plot.
#' @param obj The lcmsPlot object.
#' @param plot_config The plot configuration.
#' A list that specifies the plot function to use for each dataset type.
#' @return The patchwork plot object.
plot_multiple_datasets <- function(datasets, obj, plot_config) {
  datasets_plots <- lapply(names(datasets), function (ds_name) {
    plot_config[[ds_name]](
      datasets = datasets,
      supporting_datasets = list(detected_peaks = obj@data@detected_peaks),
      options = obj@options
    )
  })
  names(datasets_plots) <- toupper(substr(names(datasets), 1, 1))

  plt <- patchwork::wrap_plots(datasets_plots, ncol = 1, design = obj@options$layout$design)
  return(plt)
}

#' Plot multiple faceted dataset types.
#'
#' @param datasets A list of data frames each corresponding to a dataset to plot.
#' @param obj The lcmsPlot object.
#' @param plot_config The plot configuration.
#' A list that specifies the plot function to use for each dataset type.
#' @return The patchwork plot object.
plot_multiple_faceted_datasets <- function(datasets, obj, plot_config) {
  facets_options <- obj@options$facets
  facets <- facets_options$facets

  metadata <- left_join(
    obj@data@additional_metadata,
    obj@data@metadata,
    by = c("metadata_index" = "sample_index")
  )

  all_plots <- metadata %>%
    group_by(!!!syms(facets)) %>%
    group_map(~ {
      metadata_df <- .y

      datasets_subset <- lapply(datasets, function (ds) {
        ds %>%
          filter(across(all_of(facets), ~ . == metadata_df[[cur_column()]]))
      })

      datasets_plots <- lapply(names(datasets_subset), function (ds_name) {
        plot_config[[ds_name]](
          datasets = datasets_subset,
          supporting_datasets = list(detected_peaks = obj@data@detected_peaks),
          options = obj@options
        )
      })

      plot_title <- facets %>%
        lapply(function (x) metadata_df[[x]]) %>%
        unlist() %>%
        paste0(collapse = "\n")

      datasets_plots[[1]] <- datasets_plots[[1]] +
        ggtitle(plot_title) +
        theme(plot.title = element_text(size = 9, hjust = 0.5))

      combined_plot <- patchwork::wrap_plots(
        datasets_plots,
        ncol = 1,
        design = obj@options$layout$design
      )

      combined_plot <- combined_plot + patchwork::plot_layout(axes = "collect")

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

  plt <- patchwork::wrap_plots(
    plot_grid,
    byrow = TRUE,
    ncol = ncol,
    nrow = nrow
  )

  plt <- plt + patchwork::plot_layout(axes = "collect", guides = "collect")

  return(plt)
}

plot_multiple_gridded_datasets <- function(datasets, obj, plot_config) {
  stop("Currently not supported")
}
