#' Plot a single dataset type
#'
#' `plot_single_dataset()` allows the plotting of a single dataset type
#' (e.g., only chromatograms).
#'
#' @param datasets A named list of data frames, where each element represents
#' a dataset to be plotted. The names must match the supported dataset types
#' defined in `plot_config`.
#' @param obj An instance of class `lcmsPlotClass`.
#' @param plot_config A `list` that defines the plot configuration.
#' This list should use the supported dataset types as keys
#' (e.g., `chromatograms`), with each value providing the corresponding
#' function used to plot that dataset type.
#' @return A `patchwork` object combining all generated plots. In this
#' case the `patchwork` object contains a single `ggplot2` object.
#' @keywords internal
plot_single_dataset <- function(datasets, obj, plot_config) {
    dataset_name <- names(datasets)

    plt <- plot_config[[dataset_name]](
        datasets = datasets,
        supporting_datasets = list(
                detected_peaks = obj@data@detected_peaks,
                purity_scores = obj@data@purity_scores),
        options = obj@options,
        single = TRUE
    )

    plt <- patchwork::wrap_plots(plt)
    return(plt)
}

#' Plot multiple dataset types
#'
#' `plot_multiple_datasets()` iterates over the provided datasets, dispatches
#' each one to its corresponding plotting function as defined in
#' `plot_config`, and combines the resulting plots into a single patchwork
#' object.
#'
#' @param datasets A named `list` of data frames, where each element represents
#' a dataset to be plotted. The names must match the supported dataset types
#' defined in `plot_config`.
#' @param obj An instance of class `lcmsPlotClass`.
#' @param plot_config A named `list` defining the plotting functions to use for
#' each dataset type. List names correspond to dataset types, and values are
#' functions that return ggplot objects.
#' @return A `patchwork` object combining all generated plots.
#' @keywords internal
plot_multiple_datasets <- function(datasets, obj, plot_config) {
    datasets_plots <- lapply(names(datasets), function (ds_name) {
        plot_config[[ds_name]](
            datasets = datasets,
            supporting_datasets = list(
                detected_peaks = obj@data@detected_peaks),
            options = obj@options
        )
    })
    names(datasets_plots) <- toupper(substr(names(datasets), 1, 1))

    plt <- patchwork::wrap_plots(
        datasets_plots,
        ncol = 1,
        design = obj@options$layout$design)
    return(plt)
}

#' Plot multiple faceted dataset types
#'
#' `plot_multiple_faceted_datasets()` generates plots for multiple dataset types
#' and arranges them into a faceted grid. Datasets are split according to the
#' faceting variables defined in `obj@options$facets`, with each facet panel
#' containing a vertically stacked set of dataset-specific plots.
#'
#' @param datasets A named `list` of data frames, where each element represents
#' a dataset to be plotted. The names must match the supported dataset types
#' defined in `plot_config`.
#' @param obj An instance of class `lcmsPlotClass`.
#' @param plot_config A named `list` defining the plotting functions to use for
#' each dataset type. List names correspond to dataset types, and values are
#' functions that return ggplot objects.
#' @return A `patchwork` object combining all generated plots.
#' @keywords internal
plot_multiple_faceted_datasets <- function(datasets, obj, plot_config) {
    facets_options <- obj@options$facets
    facets <- facets_options$facets

    if (nrow(obj@data@feature_metadata) == 0) {
        metadata <- obj@data@metadata
    } else {
        metadata <- left_join(
            obj@data@feature_metadata,
            obj@data@metadata,
            by = c("metadata_index" = "sample_index")
        )
    }

    all_plots <- metadata |>
        group_by(!!!syms(facets)) |>
        group_map(~ {
            metadata_df <- .y

            datasets_subset <- lapply(datasets, function (ds) {
                # This keeps only rows in 'ds' that match the values
                # in 'metadata_df' for the columns listed in 'facets'
                inner_join(ds, metadata_df, by = facets)
            })

            datasets_plots <- lapply(names(datasets_subset), function (dsn) {
                plot_config[[dsn]](
                    datasets = datasets_subset,
                    supporting_datasets = list(
                        detected_peaks = obj@data@detected_peaks),
                    options = obj@options
                )
            })

            plot_title <- facets |>
                lapply(function (x) metadata_df[[x]]) |>
                unlist() |>
                paste0(collapse = "\n")

            datasets_plots[[1]] <- datasets_plots[[1]] +
                ggtitle(plot_title) +
                theme(plot.title = element_text(size = 9, hjust = 0.5))

            combined_plot <- patchwork::wrap_plots(
                datasets_plots,
                ncol = 1,
                design = obj@options$layout$design
            )

            combined_plot <- combined_plot +
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
