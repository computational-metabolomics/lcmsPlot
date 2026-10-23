#' Create an `lcmsPlotClass` object
#'
#' The `lcmsPlotClass` class allows a unified approach for  the management
#' of LC-MS data for the purpose of visualisation. It includes the options
#' for customising the plot, the LC-MS data, and the underlying plot object.
#' The `lcmsPlot` function is the main entry point and the preferred approach
#' to creating `lcmsPlotClass` objects.
#'
#' @param dataset The data to plot: an `XcmsExperiment`, `MsExperiment`,
#' `XCMSnExp`, `MChromatograms`, `XChromatograms`, `XChromatogram`,
#' `XcmsRawList`, `purityA`, a data source (see `MZmineFeatureListsSource()`,
#' `MsDialPeaksSource()`, `CompoundDiscovererNodeSource()`,
#' `LipidSearchSource()`), a Compound Discoverer `.cdResult` path, or a
#' `character` vector of raw file paths (`.mzML`, `.mzXML`, `.CDF`, `.raw`).
#' @param sample_id_column A `character` value indicating which column
#' should be used as the sample ID. By default it is `"sample_id"`.
#' @param metadata A `data.frame` containing the samples metadata
#' in case it is not provided in the dataset object.
#' @param BPPARAM A `BiocParallelParam` object for enabling parallelism.
#' See \link[BiocParallel:BiocParallelParam-class]{BiocParallelParam}
#' for more information.
#' @param batch_size A `numeric` value indicating the number of
#' samples per batch. This parameter is necessary when plotting
#' multiple batches using the `iterate_plot_batches` or `next_plot` functions.
#' @return An instance of `lcmsPlotClass`. It will create the necessary
#' internal structures related to the data (`data` slot)
#' and options (`options` slot).
#' @export
#' @examples
#' raw_files <- dir(
#'    system.file("cdf", package = "faahKO"),
#'    full.names = TRUE,
#'    recursive = TRUE)[1:5]
#'
#' p <- lcmsPlot(raw_files)
lcmsPlot <- function(
    dataset,
    sample_id_column = "sample_id",
    metadata = NULL,
    batch_size = NULL,
    BPPARAM = NULL
) {
    opts <- default_options()
    opts$sample_id_column <- sample_id_column
    opts$batch_size <- batch_size
    opts$parallel_param <- BPPARAM

    new("lcmsPlotClass",
        options = opts,
        data = lcmsPlotDataContainer(
            dataset,
            sample_id_column,
            metadata),
        plot = NULL)
}

setOldClass(c("gg", "ggplot"))

#' Managing LC-MS data for visualisation
#'
#' The `lcmsPlotClass` class allows a unified approach for  the management
#' of LC-MS data for the purpose of visualisation. It includes the options
#' for customising the plot, the LC-MS data, and the underlying plot object.
#'
#' @section General information:
#' The `lcmsPlotClass` class has been designed to be the entry point for
#' all data and outputs related to the `lcmsPlot` package.
#' The class abstracts away the data handling, making it easier to use
#' `lcmsPlot` with existing data wrappers like `MsExperiment` or `XCMSnExp`.
#'
#' @section Preferred usage:
#' The `lcmsPlotClass` class can be used directly to instantiate an object,
#' however the preferred approach is to use the `lcmsPlot` function.
#'
#' @slot options A `list` to store the plot options.
#' @slot data An instance of class `lcmsPlotDataContainer`.
#' @slot history A `list` to store the applied layers to generate a plot;
#' for internal use.
#' @slot plot A `patchwork` object representing the underlying plot object.
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
        options = default_options(),
        data = NULL,
        history = list(),
        plot = NULL
    )
)

#' Apply a function to an `lcmsPlotClass` object using the infix `+` operator
#'
#' This provides a convenient infix style for applying transformations to
#' `lcmsPlotClass` objects.
#'
#' @param e1 An instance of class `lcmsPlotClass`.
#' @param e2 A function that takes an `lcmsPlotClass` object
#' and returns another.
#' @return An instance of class `lcmsPlotClass`.
#' @export
#' @examples
#' raw_files <- dir(
#'    system.file("cdf", package = "faahKO"),
#'    full.names = TRUE,
#'    recursive = TRUE)[1:5]
#'
#' p <- lcmsPlot(raw_files) +
#'   lp_chromatogram(aggregation_fun = "max") +
#'   lp_arrange(group_by = "sample_id") +
#'   lp_legend(position = "bottom") +
#'   lp_labels(legend = "Sample")
setMethod(
    f = "+",
    signature = c("lcmsPlotClass", "function"),
    definition = function(e1, e2) { e2(e1) }
)

#' Move to the next plot object (batch-mode)
#'
#' `next_plot` progresses an `lcmsPlotClass` object to the next plot in a
#' batch-processing sequence.
#' This is typically used when multiple plots are generated and inspected
#' iteratively, such as when navigating large LC–MS datasets
#' in a batched workflow. The batch size is defined in the `lcmsPlot`
#' function's argument `batch_size`.
#'
#' @param object An instance of class `lcmsPlotClass`.
#' @export
#' @return An instance of class `lcmsPlotClass`.
#' @examples
#' raw_files <- dir(
#'    system.file("cdf", package = "faahKO"),
#'    full.names = TRUE,
#'    recursive = TRUE)[1:5]
#'
#' p <- lcmsPlot(raw_files, batch_size = 2) +
#'   lp_chromatogram(features = rbind(c(
#'     mzmin = 334.9,
#'     mzmax = 335.1,
#'     rtmin = 2700,
#'     rtmax = 2900))) +
#'   lp_arrange(group_by = "sample_id") +
#'   lp_legend(position = "bottom") +
#'   lp_labels(legend = "Sample")
#'
#' p <- next_plot(p)
#' p
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
#' `iterate_plot_batches` iterates over batches of plots defined by the
#' `batch_size` parameter passed to the `lcmsPlot` constructor function.
#'
#' @param object An instance of class `lcmsPlotClass`.
#' @param iter_fn The function to apply to each item being iterated on.
#' @return \code{NULL} (called for its side effect).
#' @export
#' @examples
#' raw_files <- dir(
#'    system.file("cdf", package = "faahKO"),
#'    full.names = TRUE,
#'    recursive = TRUE)[1:5]
#'
#' p <- lcmsPlot(raw_files, batch_size = 2) +
#'   lp_chromatogram(features = rbind(c(
#'     mzmin = 334.9,
#'     mzmax = 335.1,
#'     rtmin = 2700,
#'     rtmax = 2900))) +
#'   lp_arrange(group_by = "sample_id") +
#'   lp_legend(position = "bottom") +
#'   lp_labels(legend = "Sample")
#'
#' pdf(tempfile(fileext = ".pdf"))
#' iterate_plot_batches(p, function(plot_obj) {
#'   print(plot_obj)
#' })
#' dev.off()
setGeneric(
    "iterate_plot_batches",
    function(object, iter_fn) standardGeneric("iterate_plot_batches")
)

#' @rdname iterate_plot_batches
setMethod(
    f = "iterate_plot_batches",
    signature = c("lcmsPlotClass", "function"),
    function(object, iter_fn) {
        if (is.null(object@options$batch_size)) {
            stop("iterate_plot_batches requires batch_size")
        }

        sample_ids <- sampleMetadata(object)$sample_id

        if (length(sample_ids) > object@options$batch_size) {
            split_f <- ceiling(
                seq_along(sample_ids) / object@options$batch_size)
            batches <- split(sample_ids, split_f)
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
    }
)

#' Plot the `lcmsPlotClass` object
#'
#' Display an instance of `lcmsPlotClass` class to the selected device.
#'
#' @param object An instance of class `lcmsPlotClass`.
#' @return Invisible \code{NULL}
#' @export
#' @examples
#' raw_files <- dir(
#'    system.file("cdf", package = "faahKO"),
#'    full.names = TRUE,
#'    recursive = TRUE)[1:5]
#'
#' ## Shows summary information as the plot has not been built yet
#' p <- lcmsPlot(raw_files)
#' p
#'
#' ## Shows the actual plot
#' p <- lcmsPlot(raw_files) +
#'   lp_chromatogram(aggregation_fun = "max") +
#'   lp_arrange(group_by = "sample_id") +
#'   lp_legend(position = "bottom") +
#'   lp_labels(legend = "Sample")
#'
#' p
setMethod(
    f = "show",
    signature = "lcmsPlotClass",
    function(object) {
        if (.has_data(object)) {
            if (!object@options$bypass_plot_generation) {
                object <- .render_plot(
                    object,
                    additional_datasets = .additional_datasets(object))
            }
            object@options$bypass_plot_generation <- FALSE
            print(object@plot)
        } else {
            sample_data <- sampleMetadata(object)
            cat("Object of class", class(object), "\n")
            cat(" Source type:", class(object@data@source)[1], "\n")
            cat(
                " Metadata:",
                paste0(
                    nrow(sample_data),
                    " rows, ",
                    ncol(sample_data),
                    " columns"
                ),
                "\n"
            )
            cat(" Sample ID column:", object@options$sample_id_column, "\n")
            cat(" NOTE: No data has been requested to plot.\n")
        }
    }
)

make_interface_function <- function(name, args_list, fn) {
    function(obj, record_history = TRUE) {
        if (record_history) {
            obj@history <- c(
                obj@history,
                list(list(name = name, args = args_list)))
        }

        fn(obj)
    }
}

#' Get the underlying plot object.
#'
#' @return A function that takes an `lcmsPlot` object and returns a modified
#' version with the rendered plot stored in the `plot` slot.
#' It is intended for use with the `+` operator, which incrementally layers
#' new data or visual components onto the `lcmsPlot` object.
#' @export
#' @examples
#' raw_files <- dir(
#'    system.file("cdf", package = "faahKO"),
#'    full.names = TRUE,
#'    recursive = TRUE)[1:4]
#'
#' ## Create faceted chromatogram plots with a reference RT line
#' p <- lcmsPlot(raw_files) +
#'   lp_chromatogram(features = rbind(c(
#'     mzmin = 334.9,
#'     mzmax = 335.1,
#'     rtmin = 2700,
#'     rtmax = 2900))) +
#'   lp_facets(facets = 'sample_id', ncol = 4) +
#'   lp_rt_line(intercept = 2800, line_type = 'solid', color = 'red')
#' p
#'
#' ## Extract the ggplot object and apply a theme
#' p <- p +
#'   lp_get_plot() +
#'   ggplot2::theme_bw()
#' p
lp_get_plot <- function() {
    function(obj) {
        obj <- .render_plot(
            obj,
            additional_datasets = .additional_datasets(obj))
        return(obj@plot)
    }
}

#' @rdname lcmsPlotDataContainer-accessors
#' @export
setMethod("sampleMetadata", "lcmsPlotClass", function(object, ...) {
    sampleMetadata(object@data)
})

#' @rdname lcmsPlotDataContainer-accessors
#' @export
setMethod("spectra", "lcmsPlotClass", function(object, ...) {
    spectra(object@data)
})

#' @rdname lcmsPlotDataContainer-accessors
#' @export
setMethod("chromatograms", "lcmsPlotClass", function(object, ...) {
    chromatograms(object@data)
})

#' @rdname lcmsPlotDataContainer-accessors
#' @export
setMethod("detectedPeaks", "lcmsPlotClass", function(object, ...) {
    detectedPeaks(object@data)
})

#' @rdname lcmsPlotDataContainer-accessors
#' @export
setMethod("summaries", "lcmsPlotClass", function(object, ...) {
    summaries(object@data)
})
