#' Create an lcmsPlotClass object
#'
#' @param dataset An object of type XCMSnExp, MsExperiment, or character
#' @param sample_id_column Which column should be used as the sample ID
#' @param metadata The metadata in case it's not provided in the dataset object
#' @param parallel_param The BiocParallel object for enabling parallelism
#' @returns An lcmsPlotClass object
#' @export
lcmsPlot <- function(
    dataset,
    sample_id_column = "sample_id",
    metadata = NULL,
    parallel_param = NULL,
    batch_size = NULL
) {
  opts <- default_options
  opts$sample_id_column <- sample_id_column
  opts$parallel_param <- parallel_param
  opts$batch_size <- batch_size

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
    history = "list",
    plot = "ANY"
  ),
  prototype = list(
    options = default_options,
    data = NULL,
    history = list(),
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

#' Set the plot for an lcmsPlotClass object
#'
#' @param object The lcmsPlotClass object
#' @param additional_datasets Additional datasets to include in the plotting
#' @export
setGeneric(
  "set_plot",
  function(object, additional_datasets) standardGeneric("set_plot")
)

#' @rdname set_plot
setMethod(
  f = "set_plot",
  signature = c("lcmsPlotClass", "list"),
  function(object, additional_datasets) {
    dataset_types <- c(DATASET_TYPES, names(additional_datasets))
    datasets <- lapply(dataset_types, function(dataset_name) {
      if (object@options[[dataset_name]]$show) {
        if (dataset_name %in% slotNames(object@data)) {
          data_df <- slot(object@data, dataset_name)
        } else {
          data_df <- additional_datasets[[dataset_name]]
        }

        if (nrow(data_df) == 0) {
          stop(paste0("Empty dataset ", dataset_name))
        }

        data_df <- merge_by_index(data_df, object@data@metadata, index_col = 'metadata_index')

        if (nrow(object@data@additional_metadata) > 0) {
          data_df <- merge_by_index(data_df, object@data@additional_metadata, index_col = 'additional_metadata_index')
        }

        return(data_df)
      } else {
        return(NULL)
      }
    })
    names(datasets) <- dataset_types
    datasets <- remove_null_elements(datasets)

    object@plot <- plot_data(datasets, object)

    return(object)
  }
)

#' Set the plot for an lcmsPlotClass object
#'
#' @param object The lcmsPlotClass object
#' @export
setGeneric(
  "next_plot",
  function(object) standardGeneric("next_plot")
)

#' @rdname next_plot
setMethod(
  f = "next_plot",
  signature = c("lcmsPlotClass"),
  function(object) {
    object@options$batch_index <- object@options$batch_index + 1
    for (history_item in object@history) {
      fn <- get(history_item$name, asNamespace("lcmsPlot"))
      object <- do.call(fn, history_item$args)(object)
    }
    return(object)
  }
)

#' Iterate on the batches of plots
#'
#' @param object The lcmsPlotClass object
#' @param iter_fn The function to apply to each item being iterated on
#' @export
setGeneric(
  "iterate_plot_batches",
  function(object, iter_fn) standardGeneric("iterate_plot_batches")
)

#' @rdname iterate_plot_batches
setMethod(
  f = "iterate_plot_batches",
  signature = c("lcmsPlotClass", "function"),
  function(object, iter_fn) {
    # 1. Get sample IDs 
    
    # object@options$batch_index <- object@options$batch_index + 1
    # for (history_item in object@history) {
    #   fn <- get(history_item$name, asNamespace("lcmsPlot"))
    #   object <- do.call(fn, history_item$args)(object)
    # }
    
    if (is.null(object@options$batch_size)) {
      stop("iterate_plot_batches requires batch_size")
    }
    
    # TODO: check this
    sample_ids <- object@data@metadata$sample_id # object@options$chromatograms$sample_ids
    
    if (length(sample_ids) > object@options$batch_size) {
      batches <- split(sample_ids, ceiling(seq_along(sample_ids) / object@options$batch_size))
    } else {
      batches <- list(sample_ids)
    }
    
    object@options$batch_index <- 1
    for (batch in batches) {
      for (history_item in object@history) {
        fn <- get(history_item$name, asNamespace("lcmsPlot"))
        object <- do.call(fn, history_item$args)(object, FALSE)
      }
      iter_fn(object)
      object@options$batch_index <- object@options$batch_index + 1
    }
    
    # return(object)
  }
)

#' Plot the lcmsPlotClass object
#'
#' @param object The lcmsPlotClass object
#' @export
setMethod(
  f = "show",
  signature = "lcmsPlotClass",
  function(object) {
    if (!object@options$bypass_plot_generation) {
      object <- set_plot(object, additional_datasets = list())
    }
    object@options$bypass_plot_generation <- FALSE
    print(object@plot)
  }
)

make_interface_function <- function(name, args_list, fn) {
  function(obj, record_history = TRUE) {
    if (record_history) {
      obj@history <- c(obj@history, list(list(name = name, args = args_list)))
    }
    
    fn(obj)
  }
}

#' Define the chromatograms to plot.
#'
#' @param features A character vector or a matrix of mz and rt representing features to plot.
#' @param sample_ids A character vector of sample IDs to plot.
#' @param ppm The ppm error for the chromatograms.
#' @param rt_tol The RT tolerance for the chromatograms.
#' @param highlight_peaks Whether to highlight the picked peaks.
#' @param highlight_peaks_color The color of the highlighted peaks. By default it colors by sample.
#' @param highlight_peaks_factor The factor that determines the color.
#' @param aggregation_fun In case of plotting the full RT range, which aggregation function to use for the spectra intensities.
#' @param rt_adjusted Whether to plot the RT adjusted version of the chromatograms.
#' @param rt_unit The unit to use for the RT (one of "minute" or "second").
#' @param intensity_unit The unit to use for the intensity (one of "absolute" or "relative").
#' @param fill_gaps Whether to fill gaps in RT with 0 intensity.
#' @returns A function that takes and returns a lcmsPlotClass object.
#' @export
chromatogram <- function(
  features = NULL,
  sample_ids = NULL,
  ppm = 10,
  rt_tol = 10,
  highlight_peaks = FALSE,
  highlight_peaks_color = NULL,
  highlight_peaks_factor = "sample_id",
  aggregation_fun = "max",
  rt_adjusted = FALSE,
  rt_unit = "second",
  intensity_unit = "absolute",
  fill_gaps = FALSE
) {
  make_interface_function(
    name = "chromatogram",
    args_list = as.list(environment()),
    fn = function(obj) {
      if (is.null(sample_ids)) {
        sample_ids <- obj@data@metadata$sample_id
      }
      
      if (!is.null(obj@options$batch_size) && length(sample_ids) > obj@options$batch_size) {
        batches <- split(sample_ids, ceiling(seq_along(sample_ids) / obj@options$batch_size))
        batch_sample_ids <- batches[[obj@options$batch_index]]
      } else {
        batch_sample_ids <- sample_ids
      }
      
      obj@options$chromatograms <- list(
        show = TRUE,
        features = features,
        sample_ids = batch_sample_ids,
        ppm = ppm,
        rt_tol = rt_tol,
        highlight_peaks = highlight_peaks,
        highlight_peaks_color = highlight_peaks_color,
        highlight_peaks_factor = highlight_peaks_factor,
        aggregation_fun = aggregation_fun,
        rt_adjusted = rt_adjusted,
        rt_unit = rt_unit,
        intensity_unit = intensity_unit,
        fill_gaps = fill_gaps
      )
      
      if (is.null(features)) {
        obj@data <- create_full_rt_chromatograms(obj@data, obj@options)
      } else if (is.character(features)) {
        obj@data <- create_chromatograms_from_feature_ids(obj@data, obj@options)
      } else {
        obj@data <- create_chromatograms_from_features(obj@data, obj@options)
      }
      
      return(obj)
    }
  )
}

#' Define the mass trace to plot.
#'
#' @returns A function that takes and returns a lcmsPlotClass object.
#' @export
mass_trace <- function() {
  make_interface_function(
    name = "mass_trace",
    args_list = list(),
    fn = function(obj) {
      obj@options$mass_traces$show <- TRUE
      return(obj)
    }
  )
}

#' Define the spectra to plot
#'
#' @param sample_ids The sample IDs to consider. If NULL it will use the chromatogram ones.
#' @param mode The method to choose the scan. One of: closest, closest_apex, across_peak.
#' @param ms_level The MS level to consider for the scan.
#' @param rt The RT to consider - mode=closest
#' @param scan_index The scan index to consider
#' @param interval The RT interval to consider - mode=across_peak
#' @param spectral_match_db The spectral database to match against
#' @param match_target_index The target index for the mirror plot (index from the highest scoring)
#' @returns A function that takes and returns a lcmsPlotClass object
#' @export
spectra <- function(
  sample_ids = NULL,
  mode = 'closest_apex',
  ms_level = 1,
  rt = NULL,
  scan_index = NULL,
  interval = 3,
  spectral_match_db = NULL,
  match_target_index = NULL
) {
  make_interface_function(
    name = "spectra",
    args_list = as.list(environment()),
    fn = function(obj) {
      is_standalone <- !obj@options$chromatograms$show
      
      if  (is_standalone) {
        if (is.null(sample_ids)) {
          sample_ids <- obj@data@metadata$sample_id
        }
      } else {
        sample_ids <- obj@options$chromatograms$sample_ids
      }
      
      obj@options$spectra <- list(
        show = TRUE,
        sample_ids = sample_ids,
        mode = mode,
        ms_level = ms_level,
        rt = rt,
        scan_index = scan_index,
        interval = interval,
        spectral_match_db = spectral_match_db,
        match_target_index = match_target_index
      )
      
      obj@data <- create_spectra(obj@data, obj@options)
      return(obj)
    }
  )
}

#' Define the total ion currents.
#'
#' @param sample_ids The sample IDs to select.
#' @param type The type of plot; one of "boxplot", "violin", "jitter".
#' @export
total_ion_current <- function(sample_ids = NULL, type = "boxplot") {
  function(obj) {
    if (!inherits(obj@data@data_obj, c("XCMSnExp", "MsExperiment"))) {
      stop("total_ion_current: to plot the total ion current the data object should be either of class XCMSnExp or MsExperiment.")
    }

    if (is.null(sample_ids)) {
      sample_ids <- obj@data@metadata$sample_id
    }

    obj@options$total_ion_current <- list(
      show = TRUE,
      sample_ids = sample_ids,
      type = type
    )

    obj@data <- create_total_ion_current(obj@data, obj@options)

    return(obj)
  }
}

#' Define a 2D intensity map
#'
#' @param mz_range The m/z range of the map
#' @param rt_range The RT range of the map
#' @param sample_ids The sample IDs to select
#' @param density Whether to show a density or a point-cloud plot
#' @returns A function that takes and returns a lcmsPlotClass object
#' @export
intensity_map <- function(mz_range, rt_range, sample_ids = NULL, density = FALSE) {
  function(obj) {
    if (is.null(sample_ids)) {
      sample_ids <- obj@data@metadata$sample_id
    }

    obj@options$intensity_maps <- list(
      show = TRUE,
      sample_ids = sample_ids,
      mz_range = mz_range,
      rt_range = rt_range,
      density = density
    )

    obj@data <- create_intensity_map(obj@data, obj@options)

    return(obj)
  }
}

#' Define the RT difference plot between raw and adjusted datasets
#'
#' @returns A function that takes and returns a lcmsPlotClass object
#' @export
rt_diff_plot <- function() {
  function(obj) {
    if (!inherits(obj@data@data_obj, c("XCMSnExp", "MsExperiment"))) {
      stop("rt_diff_plot: to plot the RT differences the data object should be either of class XCMSnExp or MsExperiment.")
    }

    if (!xcms_utils$has_rt_alignment_been_performed(obj@data@data_obj)) {
      stop("rt_diff_plot: RT alignment was not performed.")
    }

    obj@options$rt_diff <- list(show = TRUE)
    obj@data <- create_rt_diff(obj@data, obj@options)

    return(obj)
  }
}

#' Define the arrangement of chromatograms
#'
#' @param group_by The column to group by (in the samples metadata)
#' @returns A function that takes and returns a lcmsPlotClass object
#' @export
arrange <- function(group_by) {
  make_interface_function(
    name = "arrange",
    args_list = as.list(environment()),
    fn = function(obj) {
      obj@options$arrangement <- list(
        group_by = group_by
      )
      return(obj)
    }
  )
}

#' Define plot's faceting
#'
#' @param facets The facet factors from the sample metadata
#' @param ncol The number of columns
#' @param nrow The number of rows
#' @param free_x Allow scales to vary across x
#' @param free_y Allow scales to vary across y
#' @returns A function that takes and returns a lcmsPlotClass object
#' @export
facets <- function(facets, ncol = NULL, nrow = NULL, free_x = FALSE, free_y = FALSE) {
  make_interface_function(
    name = "facets",
    args_list = as.list(environment()),
    fn = function(obj) {
      obj@options$facets <- list(
        facets = facets,
        ncol = ncol,
        nrow = nrow,
        free_x = free_x,
        free_y = free_y
      )
      return(obj)
    }
  )
}

#' Define a gridded plot
#'
#' @param rows The factors that represent rows
#' @param cols The factors that represent columns
#' @param free_y Whether the y-axis is free for each row
#' @returns A function that takes and returns a lcmsPlotClass object
#' @export
grid <- function(rows, cols, free_y = FALSE) {
  make_interface_function(
    name = "grid",
    args_list = as.list(environment()),
    fn = function(obj) {
      obj@options$grid <- list(
        rows = rows,
        cols = cols,
        free_y = free_y
      )
      return(obj)
    }
  )
}

#' Define the labels of the plot, such as title and legend
#'
#' @param title The plot title
#' @param legend The legend's title
#' @returns A function that takes and returns a lcmsPlotClass object
#' @export
labels <- function(title = NULL, legend = NULL) {
  make_interface_function(
    name = "labels",
    args_list = as.list(environment()),
    fn = function(obj) {
      obj@options$labels <- list(
        title = title,
        legend = legend
      )
      return(obj)
    }
  )
}

#' Define the legend layout
#'
#' @param position The legend position
#' @returns A function that takes and returns a lcmsPlotClass object
#' @export
legend <- function(position = NULL)  {
  make_interface_function(
    name = "legend",
    args_list = as.list(environment()),
    fn = function(obj) {
      obj@options$legend <- list(
        position = position
      )
      return(obj)
    }
  )
}

#' Define a vertical line across a retention time value
#'
#' @param intercept The x-axis intercept
#' @param line_type The line type
#' @param color The line color
#' @export
rt_line <- function(intercept, line_type = 'dashed', color = 'black') {
  make_interface_function(
    name = "rt_line",
    args_list = as.list(environment()),
    fn = function(obj) {
      rt_line_obj <- list(
        intercept = intercept,
        line_type = line_type,
        color = color
      )
      obj@options$rt_lines <- append(obj@options$rt_lines, list(rt_line_obj))
      return(obj)
    }
  )
}

#' Define the plot layout
#'
#' @param design Specification of the location of areas in the layout (see https://patchwork.data-imaginist.com/reference/wrap_plots.html#arg-design)
#' @export
layout <- function(design = NULL) {
  make_interface_function(
    name = "layout",
    args_list = as.list(environment()),
    fn = function(obj) {
      obj@options$layout = list(
        design = design
      )
      return(obj)
    }
  )
}

#' Get the underlying plot object
#' @export
get_plot <- function() {
  function(obj) {
    obj <- set_plot(obj, additional_datasets = list())
    return(obj@plot)
  }
}
