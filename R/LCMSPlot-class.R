#' Create a lcmsPlotClass object
#'
#' @param dataset An object of type XCMSnExp or MsExperiment
#' @param sample_id_column Which column should be used as the sample ID
#' @returns A lcmsPlotClass object
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
    dataset_types <- c("chromatograms", "mass_traces")
    datasets <- lapply(dataset_types, function(dataset_name) {
      if (object@options[[dataset_name]]$show) {
        dataset <- slot(object@data, dataset_name)
        dataset <- merge_by_index(dataset, object@data@processed_data_info, index_col = 'metadata_index')
        return(dataset)
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
#' @param features A vector of feature IDs to plot (TODO: change)
#' @param sample_ids A vector of sample IDs to plot
#' @param ppm The ppm error for the chromatograms
#' @param rt_tol The RT tolerance for the chromatograms
#' @param highlight_peaks Whether to highlight the picked peaks
#' @returns A function that takes and returns a lcmsPlotClass object
#' @export
chromatogram <- function(
  features,
  sample_ids = NULL,
  ppm = 10,
  rt_tol = 10,
  highlight_peaks = FALSE
) {
  function(obj) {
    if (is.null(sample_ids)) {
      sample_ids <- obj@data@metadata$sample_id
    }

    if (is.character(features)) {
      obj@data <- create_chromatograms_from_features(
        obj@data,
        features,
        sample_ids,
        ppm,
        rt_tol)

      obj@options$chromatograms$highlight_peaks <- highlight_peaks
    } else {
      obj@data <- create_chromatograms_from_raw(
        obj@data,
        features,
        sample_ids,
        ppm,
        rt_tol)
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

#' Define the arrangement of chromatograms
#'
#' @param rows ...
#' @param cols ...
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

#' Define the arrangement of chromatograms
#'
#' @param facets ...
#' @param ncol ...
#' @param nrow ...
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
