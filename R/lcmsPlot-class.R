#' Create an lcmsPlotClass object
#'
#' @param dataset An object of type XCMSnExp, MsExperiment, or character
#' @param sample_id_column Which column should be used as the sample ID
#' @param metadata The metadata in case it's not provided in the dataset object
#' @returns An lcmsPlotClass object
#' @export
lcmsPlot <- function(dataset, sample_id_column = "sample_id", metadata = NULL) {
  opts <- default_options
  opts$sample_id_column <- sample_id_column

  new("lcmsPlotClass",
      options = opts,
      data = create_data_container_from_obj(dataset, sample_id_column, metadata),
      plot = NULL)
}

setOldClass(c("gg", "ggplot"))

#' lcmsPlotClass class
#'
#' @slot options The object options
#' @slot data The lcmsPlotDataContainer object
#' @slot plot The underlying plot object
#' @export
setClass(
  "lcmsPlotClass",
  slots = list(
    options = "list",
    data = "lcmsPlotDataContainer",
    plot = "ANY"
  ),
  prototype = list(
    options = default_options,
    data = NULL,
    plot = NULL
  )
)

#' Apply a function to the lcmsPlotClass object
#'
#' @param e1 A lcmsPlotClass object
#' @param e2 A function that takes a lcmsPlotClass object and returns another
#' @export
setMethod(
  f = "+",
  signature = c("lcmsPlotClass", "function"),
  definition = function(e1, e2) { e2(e1) }
)

#' Plot the lcmsPlotClass object
#'
#' @param object The lcmsPlotClass object
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_ribbon facet_wrap
#' @export
setMethod(
  f = "show",
  signature = "lcmsPlotClass",
  function(object) {
    # TODO: do validation on data types
    dataset_types <- c("chromatograms", "mass_traces", "spectra", "intensity_maps")
    datasets <- lapply(dataset_types, function(dataset_name) {
      if (object@options[[dataset_name]]$show) {
        data_df <- slot(object@data, dataset_name)

        if (nrow(obj@data@processed_data_info) > 0) {
          metadata_to_merge <- obj@data@processed_data_info
        } else {
          metadata_to_merge <- object@data@metadata
        }

        data_df <- merge_by_index(data_df, metadata_to_merge, index_col = 'metadata_index')
        return(data_df)
      } else {
        return(NULL)
      }
    })
    names(datasets) <- dataset_types
    datasets <- remove_null_elements(datasets)

    object@plot <- plot_data(datasets, object)

    print(object@plot)
  }
)

#' Define the chromatograms to plot
#'
#' @param features A character vector or a matrix of mz and rt representing features to plot
#' @param sample_ids A character vector of sample IDs to plot
#' @param ppm The ppm error for the chromatograms
#' @param rt_tol The RT tolerance for the chromatograms
#' @param highlight_peaks Whether to highlight the picked peaks
#' @param highlight_peaks_color The color of the highlighted peaks. By default it colors by sample.
#' @returns A function that takes and returns a lcmsPlotClass object
#' @export
chromatogram <- function(
  features = NULL,
  sample_ids = NULL,
  ppm = 10,
  rt_tol = 10,
  highlight_peaks = FALSE,
  highlight_peaks_color = NULL,
  aggregation_fun = "max"
) {
  function(obj) {
    if (is.null(sample_ids)) {
      sample_ids <- obj@data@metadata$sample_id
    }

    obj@options$sample_ids <- sample_ids

    if (is.null(features)) {
      obj@data <- create_full_rt_chromatograms(
        obj@data,
        sample_ids,
        aggregation_fun)
    } else if (is.character(features)) {
      obj@data <- create_chromatograms_from_feature_ids(
        obj@data,
        features,
        sample_ids,
        ppm,
        rt_tol)

      obj@options$chromatograms$highlight_peaks <- highlight_peaks
      obj@options$chromatograms$highlight_peaks_color <- highlight_peaks_color
    } else {
      obj@data <- create_chromatograms_from_features(
        obj@data,
        features,
        sample_ids,
        ppm,
        rt_tol)

      obj@options$chromatograms$highlight_peaks <- highlight_peaks
      obj@options$chromatograms$highlight_peaks_color <- highlight_peaks_color
    }

    return(obj)
  }
}

#' Define the mass trace to plot
#'
#' @returns A function that takes and returns a lcmsPlotClass object
#' @export
mass_trace <- function() {
  function(obj) {
    obj@options$mass_traces$show <- TRUE
    return(obj)
  }
}

#' Define the spectra to plot
#'
#' @param sample_ids The sample IDs to consider. If NULL it will use the chromatogram ones.
#' @param mode The method to choose the scan. One of: closest, closest_apex, across_peak.
#' @param ms_level The MS level to consider for the scan.
#' @param rt The RT to consider - mode=closest
#' @param interval The RT interval to consider - mode=across_peak
#' @returns A function that takes and returns a lcmsPlotClass object
#' @export
spectra <- function(sample_ids = NULL, mode = 'closest_apex', ms_level = 1, rt = NULL, interval = 3) {
  function(obj) {
    obj@options$spectra <- list(
      show = TRUE,
      mode = mode,
      ms_level = ms_level,
      rt = rt,
      interval = interval
    )

    if (!is.null(sample_ids)) {
      obj@options$sample_ids <- sample_ids
    }

    obj@data <- create_spectra(obj@data, obj@options)

    return(obj)
  }
}

#' Define a 2D intensity map
#'
#' @param mz_range The m/z range of the map
#' @param rt_range The RT range of the map
#' @param density Whether to show a density or a point-cloud plot
#' @returns A function that takes and returns a lcmsPlotClass object
#' @export
intensity_map <- function(mz_range, rt_range, density = FALSE) {
  function(obj) {
    obj@options$intensity_maps <- list(
      show = TRUE,
      mz_range = mz_range,
      rt_range = rt_range,
      density = density
    )

    obj@data <- create_intensity_map(obj@data, obj@options)

    return(obj)
  }
}

#' Define the arrangement of chromatograms
#'
#' @param group_by The column to group by (in the samples metadata)
#' @returns A function that takes and returns a lcmsPlotClass object
#' @export
arrange <- function(group_by) {
  function(obj) {
    obj@options$arrangement <- list(
      group_by = group_by
    )
    return(obj)
  }
}

#' Define plot's faceting
#'
#' @param facets The facet factors from the sample metadata
#' @param ncol The number of columns
#' @param nrow The number of rows
#' @returns A function that takes and returns a lcmsPlotClass object
#' @export
facets <- function(facets, ncol = NULL, nrow = NULL) {
  function(obj) {
    obj@options$facets <- list(
      facets = facets,
      ncol = ncol,
      nrow = nrow
    )
    return(obj)
  }
}

#' Define a gridded plot
#'
#' @param rows The factors that represent rows
#' @param cols The factors that represent columns
#' @returns A function that takes and returns a lcmsPlotClass object
#' @export
grid <- function(rows, cols) {
  function(obj) {
    obj@options$grid <- list(
      rows = rows,
      cols = cols
    )
    return(obj)
  }
}

#' Define the labels of the plot, such as title and legend
#'
#' @param title The plot title
#' @param legend The legend's title
#' @returns A function that takes and returns a lcmsPlotClass object
#' @export
labels <- function(title = NULL, legend = NULL) {
  function(obj) {
    obj@options$labels <- list(
      title = title,
      legend = legend
    )
    return(obj)
  }
}

#' Define the legend layout
#'
#' @param position The legend position
#' @returns A function that takes and returns a lcmsPlotClass object
#' @export
legend <- function(position = NULL)  {
  function(obj) {
    obj@options$legend <- list(
      position = position
    )
    return(obj)
  }
}

#' Define a vertical line across a retention time value
#'
#' @param intercept The x-axis intercept
#' @param line_type The line type
#' @param color The line color
#' @export
rt_line <- function(intercept, line_type = 'dashed', color = 'black') {
  function(obj) {
    rt_line_obj <- list(
      intercept = intercept,
      line_type = line_type,
      color = color
    )
    obj@options$rt_lines <- append(obj@options$rt_lines, list(rt_line_obj))
    return(obj)
  }
}

#' Define the plot layout
#'
#' @param design Specification of the location of areas in the layout (see https://patchwork.data-imaginist.com/reference/wrap_plots.html#arg-design)
#' @export
layout <- function(design = NULL) {
  function(obj) {
    obj@options$layout = list(
      design = design
    )
    return(obj)
  }
}
