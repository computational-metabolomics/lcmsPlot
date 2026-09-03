#' Create an `lcmsPlotClass` object
#'
#' The `lcmsPlotClass` class allows a unified approach for  the management
#' of LC-MS data for the purpose of visualisation. It includes the options
#' for customising the plot, the LC-MS data, and the underlying plot object.
#' The `lcmsPlot` function is the main entry point and the preferred approach
#' to creating `lcmsPlotClass` objects.
#'
#' @param dataset An object of type `XCMSnExp`, `MsExperiment`,
#' `MZmineSource`, or `character`.
#' If a `character` vector is supplied, it will be interpreted as a list
#' of mzML paths.
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
        data = create_data_container_from_obj(
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

.render_plot <- function(object, additional_datasets) {
    dataset_types <- c(DATASET_TYPES, names(additional_datasets))
    datasets <- lapply(dataset_types, function(dataset_name) {
        if (object@options[[dataset_name]]$show) {
            if (dataset_name %in% slotNames(object@data)) {
                data_df <- slot(object@data, dataset_name)
            } else {
                data_df <- additional_datasets[[dataset_name]]
            }

            if (nrow(data_df) > 0) {
                data_df <- merge_by_index(
                    data_df,
                    object@data@metadata,
                    index_col = 'metadata_index'
                )

                if (nrow(object@data@feature_metadata) > 0) {
                    data_df <- data_df |>
                        left_join(
                            object@data@feature_metadata,
                            by = "feature_metadata_id"
                        )
                }

                return(data_df)
            } else {
                warning("Empty dataset ", dataset_name)
                return(NULL)
            }
        } else {
            return(NULL)
        }
    })
    names(datasets) <- dataset_types
    datasets <- remove_null_elements(datasets)

    object@plot <- plot_data(datasets, object)

    return(object)
}

.has_data <- function(object) {
    data_slots <- c(
        "chromatograms",
        "mass_traces",
        "spectra",
        "total_ion_current",
        "intensity_maps",
        "rt_diff",
        "peak_density",
        "peak_count_image",
        "purity_scores"
    )

    has_slot_data <- any(vapply(
        data_slots,
        function(s) nrow(slot(object@data, s)),
        integer(1)
    ) > 0)

    # A standalone peak-boundary panel is drawn from @detected_peaks, which is
    # not a dataset slot of its own. Detected peaks alone are not a plot
    # request, so the option has to be on as well.
    has_slot_data ||
        (isTRUE(object@options$chrom_peak_rects$show) &&
            nrow(object@data@detected_peaks) > 0)
}

# Datasets that are derived from a slot at render time rather than living in one
# of their own. `.render_plot()` appends these to DATASET_TYPES, so each name
# needs a matching options entry carrying a logical `show`.
.additional_datasets <- function(object) {
    additional <- list()

    if (nrow(object@data@purity_scores) > 0) {
        if (isTRUE(object@options$purity_timeline$show)) {
            additional$purity_timeline <- object@data@purity_scores
        }
        if (isTRUE(object@options$purity_distribution$show)) {
            additional$purity_distribution <- object@data@purity_scores
        }
    }

    # The peak rectangles own a panel only when there is no host layer to
    # decorate. Deriving that here rather than at `+` time keeps it independent
    # of the order the layers were added in.
    has_host <- isTRUE(object@options$intensity_maps$show) ||
        isTRUE(object@options$mass_traces$show)

    if (isTRUE(object@options$chrom_peak_rects$show) && !has_host &&
        nrow(object@data@detected_peaks) > 0) {
        # Project down to the peak geometry plus the join keys. @detected_peaks
        # already carries the metadata columns, and .render_plot() re-joins
        # metadata unconditionally, so anything else here would come back
        # suffixed .x/.y and break faceting and the sample filter.
        additional$chrom_peak_rects <- object@data@detected_peaks |>
            dplyr::transmute(
                mz = .data$mz,
                mzmin = .data$mzmin,
                mzmax = .data$mzmax,
                rt = .data$rt,
                rtmin = .data$rtmin,
                rtmax = .data$rtmax,
                metadata_index = .data$sample_index,
                feature_metadata_id = NA_real_
            )
    }

    additional
}


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

        # TODO: needs to be reviewed
        sample_ids <- object@data@metadata$sample_id

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
            cat("Object of class", class(object), "\n")
            cat(" Data object type:", class(object@data@data_obj), "\n")
            cat(
                " Metadata:",
                paste0(
                    nrow(object@data@metadata),
                    " rows, ",
                    ncol(object@data@metadata),
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

#' Define the chromatograms to plot
#'
#' The `lp_chromatogram` function allows the generation of different types of
#' chromatograms.
#'
#' @section Summary chromatograms:
#' In this type of chromatogram, the intensities of the spectra from each scan
#' in an LC–MS dataset are aggregated into a single value per scan.
#' To create such chromatograms do not specify the `features` parameter as
#' that will create the chromatograms for the selected features.
#' In this context, the main parameter is `aggregation_fun`, which can take
#' `max` (base peak chromatogram, BPC), `sum` (total ion current, TIC), or
#' `mean` (averaged ion chromatogram).
#'
#' @section Feature chromatograms:
#' A feature is a combination of retention time (RT) and m/z. Feature
#' chromatograms can be created by specifiying the `features` parameter.
#'
#' @param features Specifies which features to generates the chromatogram for.
#' This can be either:
#' a `matrix` with columns `mz` and `rt` (optional);
#' a `matrix` with columns `mzmin`, `mzmax`,
#' `rtmin` (optional), `rtmax` (optional);
#' a `data.frame` with columns `sample_id`, `mz` and `rt` (optional);
#' a `data.frame` with columns `sample_id`, `mzmin`, `mzmax`,
#' `rtmin` (optional), `rtmax` (optional);
#' a `character` vector representing the grouped peaks (feature) names
#' as returned by `xcms::groupnames` - requires the data to be
#' an `XCMSnExp` or `MsExperiment` object with grouped peaks.
#' @param sample_ids A `character` vector specifying the sample IDs
#' to include in the plot. If `NULL`, the function uses the sample IDs
#' specified in the `lcmsPlot` object.
#' @param ppm A `numeric` value specifying the mass accuracy (in ppm) used
#' when generating chromatograms. Ignored when the `features` parameter
#' specifies both `mzmin` and `mzmax`.
#' @param rt_tol A `numeric` value specifying the RT tolerance used
#' when generating chromatograms. Ignored when the `features` parameter
#' specifies both `rtmin` and `rtmax`.
#' @param line_type A `character` value specifying the line type (from ggplot2).
#' One of: "solid", "dashed", "dotted", "dotdash", "longdash", "twodash".
#' @param highlight_peaks A `logical` value indicating whether to highlight
#' the detected peaks; the input data must be an `XCMSnExp` or `MsExperiment`
#' object.
#' @param highlight_peaks_color A `character` value indicating the color of the
#' highlighted peaks.
#' @param highlight_peaks_mode A `character` value indicating how the peaks
#' should be highlighted. One of `"polygon"` (fills the area under the curve),
#' `"rectangle"` (draws a bounding box from rtmin to rtmax up to the peak apex),
#' or `"point"` (marks the peak apex with a point). Defaults to `"polygon"`.
#' @param highlight_peaks_factor A `character` value indicating the factor from
#' the metadata that determines the color. By default it colors by `sample_id`.
#' @param aggregation_fun A `character` value indicating which aggregation
#' function to use for the spectra intensities; one of `max` (base peak
#' chromatogram), `sum` (total ion current), or `mean` (averaged ion
#' chromatogram). Only applicable to summary chromatograms.
#' @param rt_type A `chracter` value indicating what type of RT to use for the
#' chromatograms. One of `uncorrected` (default), `corrected`, or `both`;
#' the input data must be an `XCMSnExp` or `MsExperiment` object.
#' If `both` is chosen, this will give access to a metadata column called
#' `rt_adjusted` that can be used to differentiate the two RT types (e.g.,
#' through faceting).
#' @param rt_unit A `character` value indicating the unit to use for
#' the RT axis; one of `"minute"` or `"second"`.
#' @param intensity_unit A `character` value indicating the unit to use for
#' the intensity axis; one of `"absolute"` or `"relative"`.
#' @param fill_gaps A `logical` value indicating whether to fill gaps
#' in RT with 0 intensity.
#' @param na.rm A `logical` value. When `TRUE`, data points whose intensity
#' is `NA` are removed before plotting. Defaults to `FALSE`.
#' @param highlight_apices A `logical` value indicating whether to
#' highlight apices with the corresponding RT values in a chromatogram.
#' @param stacked A `numeric` value in `[0, 1]`. When non-zero, each series in a
#' panel is offset up the y axis in m/z order by that fraction of the intensity
#' range, so co-eluting traces stop occluding each other.
#' Defaults to `0` (no stacking).
#' Which series get offset follows `lp_arrange(group_by = )`; without it,
#' samples are stacked.
#' @param transform A `function` applied to the intensities before plotting,
#' such as `log10`, to compress the dynamic range. It is applied to the peak
#' highlight geometry as well, so shaded peaks stay attached to their traces.
#' Defaults to `identity`. Note that `log10` maps the exact zeros injected by
#' `fill_gaps = TRUE` to `-Inf`.
#' @return This function returns another function that takes an `lcmsPlot`
#' object and produces a modified version containing the generated chromatograms
#' in its `data` slot. It is designed to be used with the `+` operator,
#' which serves as a layering mechanism. Each use of `+` incrementally enriches
#' the `lcmsPlot` object by adding new data or visual components.
#' @export
#' @examples
#' raw_files <- dir(
#'    system.file("cdf", package = "faahKO"),
#'    full.names = TRUE,
#'    recursive = TRUE)[1:5]
#'
#' p <- lcmsPlot(raw_files) +
#'   lp_chromatogram(aggregation_fun = "max") +
#'   lp_arrange(group_by = "sample_id")
#'
#' p
lp_chromatogram <- function(
    features = NULL,
    sample_ids = NULL,
    ppm = 10,
    rt_tol = 10,
    line_type = "solid",
    highlight_peaks = FALSE,
    highlight_peaks_color = NULL,
    highlight_peaks_mode = "polygon",
    highlight_peaks_factor = "sample_id",
    aggregation_fun = "max",
    rt_type = "uncorrected",
    rt_unit = "second",
    intensity_unit = "absolute",
    fill_gaps = FALSE,
    na.rm = FALSE,
    highlight_apices = list(column = NULL, top_n = NULL),
    stacked = 0,
    transform = identity
) {
    aggregation_fun <- match.arg(aggregation_fun, c("max", "sum", "mean"))

    if (!is.function(transform)) {
        stop("lp_chromatogram: 'transform' must be a function.")
    }
    if (!is.numeric(stacked) || length(stacked) != 1L) {
        stop("lp_chromatogram: 'stacked' must be a single numeric value.")
    }

    make_interface_function(
        name = "lp_chromatogram",
        args_list = as.list(environment()),
        fn = function(obj) {
            if (is.null(sample_ids)) {
                sample_ids <- obj@data@metadata$sample_id
            }

            if (!is.null(obj@options$batch_size) &&
                length(sample_ids) > obj@options$batch_size) {
                split_f <- ceiling(
                    seq_along(sample_ids) / obj@options$batch_size)
                batches <- split(sample_ids, split_f)
                batch_sample_ids <- batches[[obj@options$batch_index]]
            } else {
                batch_sample_ids <- sample_ids
            }

            obj@options$chromatograms <- list(
                show = TRUE,
                features = features,
                sample_ids = batch_sample_ids,
                ppm = ppm,
                line_type = line_type,
                rt_tol = rt_tol,
                highlight_peaks = highlight_peaks,
                highlight_peaks_color = highlight_peaks_color,
                highlight_peaks_mode = highlight_peaks_mode,
                highlight_peaks_factor = highlight_peaks_factor,
                aggregation_fun = aggregation_fun,
                rt_type = rt_type,
                rt_unit = rt_unit,
                intensity_unit = intensity_unit,
                fill_gaps = fill_gaps,
                highlight_apices = highlight_apices,
                stacked = stacked,
                transform = transform
            )

            result <- create_chromatograms(
                obj@data@data_obj,
                obj@data@metadata,
                obj@options,
                features)

            obj@data@chromatograms <- if (na.rm) {
                result$chromatograms[!is.na(result$chromatograms$intensity), ]
            } else {
                result$chromatograms
            }
            obj@data@mass_traces <- result$mass_traces
            obj@data@feature_metadata <- result$feature_metadata
            obj@data@detected_peaks <- result$detected_peaks
            validObject(obj@data)

            return(obj)
        }
    )
}

#' Define the mass trace to plot
#'
#' The `lp_mass_trace` function enables the generation of mass traces,
#' which are graphical representations commonly used in mass spectrometry
#' data analysis. A mass trace plots individual data points defined by their
#' retention time and corresponding mass-to-charge ratio (m/z), making it easier
#' to visualise how specific ions behave over the course of a
#' chromatographic run.
#'
#' @return This function returns another function that takes an `lcmsPlot`
#' object and produces a modified version containing the generated mass traces
#' in its `data` slot. It is designed to be used with the `+` operator,
#' which serves as a layering mechanism. Each use of `+` incrementally enriches
#' the `lcmsPlot` object by adding new data or visual components.
#' @export
#' @examples
#' raw_files <- dir(
#'    system.file("cdf", package = "faahKO"),
#'    full.names = TRUE,
#'    recursive = TRUE)[1:5]
#'
#' ## Create chromatograms of a specific feature
#' p <- lcmsPlot(raw_files) +
#'   lp_chromatogram(features = rbind(c(
#'     mzmin = 334.9,
#'     mzmax = 335.1,
#'     rtmin = 2700,
#'     rtmax = 2900))) +
#'   lp_arrange(group_by = "sample_id")
#'
#' ## Add mass traces
#' p <- p + lp_mass_trace()
#'
#' p
lp_mass_trace <- function() {
    make_interface_function(
        name = "lp_mass_trace",
        args_list = list(),
        fn = function(obj) {
            obj@options$mass_traces$show <- TRUE
            return(obj)
        }
    )
}

#' Define the spectra to plot
#'
#' The `lp_spectra` function enables the generation of spectra, which are
#' graphical representations of ions detected at each mass-to-charge ratio (m/z)
#' with their corresponding absolute or relative intensities.
#'
#' @section Spectra associated with chromatograms:
#' A spectrum is obtained from a scan at a specific retention time (RT).
#' Therefore, when plotting a chromatogram together with its associated spectra,
#' it is common to mark the RT with a vertical line on the chromatogram
#' to indicate where the spectra were acquired. See the example below on
#' how to generate these types of spectra.
#'
#' @section Standalone spectra:
#' Standalone spectra can also be generated, provided no chromatograms
#' are present (i.e., lp_chromatogram has not been used).
#'
#' @param sample_ids A `character` vector specifying the sample IDs
#' to include in the plot. If `NULL`, the function uses the sample IDs
#' specified in the `lcmsPlot` object or the `lp_chromatogram` function.
#' @param mode The method to choose the scan from which to extract the spectra.
#' One of: `closest`, the closest scan to the specified RT - `rt` parameter);
#' `closest_apex`, the closest scan to a detected peak;
#' `across_peak`, selects scans across a detected peak at a certain interval
#' specified in the `interval` parameter.
#' `mode` is not applicable to standalone spectra.
#' @param ms_level The MS level to consider for the scan.
#' @param rt When `mode = "closest"`, the RT to consider.
#' @param scan_index The exact scan index to consider for extracting a spectrum.
#' `scan_index` and `mode` are mutually exclusive.
#' @param interval When `mode = "across_peak."` The RT interval to consider.
#' @param spectral_match_db The database containing reference spectra used
#' for matching and comparison with the input spectra.
#' @param match_target_index The index, ranked by descending match score,
#' identifying which reference spectrum to display in the mirror plot.
#' @param peak_label_size A `numeric` value controlling the font size of
#' m/z labels annotated on spectral peaks. Defaults to `3`.
#' @param intensity_breaks_by A `numeric` value specifying the step size
#' (in percent) between y-axis intensity breaks. Defaults to `20`.
#' @param mz_breaks_n A `numeric` value specifying the approximate number of
#' m/z axis breaks, passed to `scales::pretty_breaks(n = ...)`.
#' Defaults to `6`.
#' @param auto_facet A `logical` value. When `TRUE` (default), a facet is
#' automatically added to separate spectra from different samples or scans.
#' Set to `FALSE` to suppress automatic faceting.
#' @return This function returns another function that takes an `lcmsPlot`
#' object and produces a modified version containing the generated spectra
#' in its `data` slot. It is designed to be used with the `+` operator,
#' which serves as a layering mechanism. Each use of `+` incrementally enriches
#' the `lcmsPlot` object by adding new data or visual components.
#' @export
#' @examples
#' raw_files <- dir(
#'    system.file("cdf", package = "faahKO"),
#'    full.names = TRUE,
#'    recursive = TRUE)[1]
#'
#' p <- lcmsPlot(raw_files) +
#'   lp_chromatogram(features = rbind(c(
#'     mzmin = 334.9,
#'     mzmax = 335.1,
#'     rtmin = 2700,
#'     rtmax = 2900))) +
#'   lp_spectra(mode = "closest", rt = 2785)
#' p
lp_spectra <- function(
    sample_ids = NULL,
    mode = 'closest_apex',
    ms_level = 1,
    rt = NULL,
    scan_index = NULL,
    interval = 3,
    spectral_match_db = NULL,
    match_target_index = NULL,
    peak_label_size = 3,
    intensity_breaks_by = 20,
    mz_breaks_n = 6,
    auto_facet = TRUE
) {
    mode <- match.arg(mode, c("closest_apex", "closest", "across_peak"))

    make_interface_function(
        name = "lp_spectra",
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
                match_target_index = match_target_index,
                peak_label_size = peak_label_size,
                intensity_breaks_by = intensity_breaks_by,
                mz_breaks_n = mz_breaks_n,
                auto_facet = auto_facet
            )

            obj@data <- create_spectra(obj@data, obj@options)
            return(obj)
        }
    )
}

#' Plot the peak density for one or more m/z features
#'
#' The `lp_peak_density` function produces a peak density plot.
#' For each supplied feature (m/z bin):
#' - The x axis shows retention time.
#' - The y axis shows sample indices (1 to n), positioned within the density
#'   range.
#' - Detected peaks are plotted as coloured points at each sample's position.
#' - The kernel density estimate of peak apex RTs is drawn as a line.
#' - When `min_fraction` is supplied, the density-descent grouping algorithm
#'   is simulated and feature groups that pass the threshold are highlighted
#'   with semi-transparent rectangles.
#'
#' @param features A `matrix` or `data.frame` with columns `mzmin` and `mzmax`
#' (required) and `rtmin`, `rtmax` (optional). Each row defines one m/z bin.
#' When omitted, the features are taken from `lp_chromatogram` if it has
#' already been called on the same object. For `XChromatograms` and
#' `XChromatogram` objects it is redundant: each extracted ion chromatogram
#' already carries its own m/z and retention-time window, which is used instead.
#' @param bw A `numeric` value specifying the kernel density bandwidth in
#' seconds. Passed to [stats::density()]. Defaults to `30`.
#' @param min_fraction A `numeric` value in `[0, 1]`. When supplied, simulates
#' the `PeakDensityParam` grouping algorithm and draws semi-transparent
#' rectangles for feature groups where at least this fraction of samples (per
#' group) contain a peak. Defaults to `NULL` (no simulation).
#' @param min_samples An `integer` specifying the minimum absolute number of
#' samples per group required to define a feature group. Only used when
#' `min_fraction` is set. Defaults to `1L`.
#' @param sample_groups A vector of length equal to the number of samples
#' assigning each sample to a group (as in `PeakDensityParam`). Defaults to
#' `NULL`, which treats all samples as a single group.
#' @param max_features An `integer` specifying the maximum number of feature
#' group rectangles to draw per m/z bin. Defaults to `50L`.
#' @param rt_unit A `character` value indicating the unit for the RT axis;
#' one of `"second"` (default) or `"minute"`.
#' @param simulate A `logical` value mirroring `xcms::plotChromPeakDensity()`.
#' `TRUE` descends the density curve to derive the feature groups the supplied
#' parameters *would* produce, using `min_fraction` (defaulting to `0.5` when
#' unset). `FALSE` instead draws the feature definitions the object already
#' stores, ignoring `min_fraction` and `min_samples`. Defaults to `NULL`, which
#' simulates whenever `min_fraction` is given.
#' @return This function returns another function that takes an `lcmsPlot`
#' object and produces a modified version containing the generated peak density
#' data in its `data` slot. It is designed to be used with the `+` operator,
#' which serves as a layering mechanism.
#' @export
#' @examples
#' data_obj <- get_XCMSnExp_object_example(
#'   indices = 1:3,
#'   should_group_peaks = TRUE)
#' p <- lcmsPlot(data_obj, sample_id_column = "sample_name") +
#'   lp_peak_density(
#'     features = rbind(c(mzmin = 334.9, mzmax = 335.1,
#'                        rtmin = 2700,  rtmax = 2900)),
#'     bw = 30,
#'     min_fraction = 0.5)
#' p
lp_peak_density <- function(
    features = NULL,
    bw = 30,
    min_fraction = NULL,
    min_samples = 1L,
    sample_groups = NULL,
    max_features = 50L,
    rt_unit = "second",
    simulate = NULL
) {
    make_interface_function(
        name = "lp_peak_density",
        args_list = as.list(environment()),
        fn = function(obj) {
            if (is.null(get_detected_peaks(obj@data@data_obj))) {
                stop(
                    "lp_peak_density: the data object reports no ",
                    "chromatographic peaks. Provide an XCMSnExp, ",
                    "MsExperiment, XChromatograms, or XChromatogram object ",
                    "with detected peaks."
                )
            }

            if (is.null(features)) {
                features <- obj@options$chromatograms$features
            }
            if (is.null(features)) {
                # Chromatogram objects already define one m/z window per EIC,
                # which makes `features` redundant for them.
                features <- get_feature_windows(obj@data@data_obj)
            }
            if (is.null(features)) {
                stop(
                    "lp_peak_density: 'features' is missing and no ",
                    "chromatogram features are set. Either provide ",
                    "'features' or call lp_chromatogram() first."
                )
            }

            features <- as.data.frame(features)
            if (!all(c("mzmin", "mzmax") %in% names(features))) {
                stop("lp_peak_density: 'features' must have columns 'mzmin' and 'mzmax'.")
            }
            if (!"rtmin" %in% names(features)) features$rtmin <- NA_real_
            if (!"rtmax" %in% names(features)) features$rtmax <- NA_real_

            obj@options$peak_density <- list(
                show = TRUE,
                features = features,
                bw = bw,
                min_fraction = min_fraction,
                min_samples = min_samples,
                sample_groups = sample_groups,
                max_features = max_features,
                rt_unit = rt_unit,
                simulate = simulate
            )

            obj@data <- create_peak_density(obj@data, obj@options)

            return(obj)
        }
    )
}

#' Plot chromatographic peak counts per retention-time bin
#'
#' `lp_peak_count_image()` produces a common quality-control view:
#' retention-time bins on the x axis, samples on
#' the y axis, and fill showing how many chromatographic peaks each sample
#' yielded in each bin.
#'
#' Samples are ordered by injection order, and bins
#' containing no peaks are drawn as zero rather than dropped.
#'
#' @param bin_size A `numeric` value giving the retention-time bin width in
#' seconds. Defaults to `30`.
#' @param log A `logical` value. When `TRUE`, counts are shown on a `log2`
#' scale, which is the usual way to read the plot when counts are skewed. Empty
#' bins are left blank rather than being drawn as `-Inf`.
#' @param sample_ids A `character` vector of sample IDs to include.
#' `NULL` uses all samples.
#' @param rt_range An optional length-2 `numeric` limiting the binned
#' retention-time range. When `NULL` (default) the full acquisition range of the
#' data object is used, so a short run shows as empty bins at the edge.
#' @param fill_scale A ggplot2 scale object to use for the fill aesthetic. When
#' `NULL` (default), uses `scale_fill_viridis_c`.
#' @return This function returns another function that takes an `lcmsPlot`
#' object and produces a modified version containing the generated peak counts
#' in its `data` slot. It is designed to be used with the `+` operator,
#' which serves as a layering mechanism.
#' @export
#' @examples
#' data_obj <- get_XCMSnExp_object_example(
#'   indices = 1:3,
#'   should_group_peaks = TRUE)
#'
#' p <- lcmsPlot(data_obj, sample_id_column = "sample_name") +
#'   lp_peak_count_image(bin_size = 30)
#' p
lp_peak_count_image <- function(
    bin_size = 30,
    log = FALSE,
    sample_ids = NULL,
    rt_range = NULL,
    fill_scale = NULL
) {
    if (!is.numeric(bin_size) || length(bin_size) != 1L || bin_size <= 0) {
        stop("lp_peak_count_image: 'bin_size' must be a positive number.")
    }

    make_interface_function(
        name = "lp_peak_count_image",
        args_list = as.list(environment()),
        fn = function(obj) {
            if (is.null(get_detected_peaks(obj@data@data_obj))) {
                stop(
                    "lp_peak_count_image: the data object reports no ",
                    "chromatographic peaks. Provide an XCMSnExp or ",
                    "MsExperiment object with detected peaks."
                )
            }

            obj@options$peak_count_image <- list(
                show = TRUE,
                sample_ids = sample_ids,
                bin_size = bin_size,
                log = log,
                rt_range = rt_range,
                fill_scale = fill_scale
            )

            obj@data <- create_peak_count_image(obj@data, obj@options)

            return(obj)
        }
    )
}

#' Define the total ion current (TIC)
#'
#' The `lp_total_ion_current` generates summary data for the
#' total ion current (TIC) of the selected samples. It works both with
#' XCMS objects (`XCMSnExp`, `MsExperiment`) and with raw files passed
#' directly to `lcmsPlot()` (a `character` vector of `.mzML`, `.mzXML`,
#' `.CDF`, or `.raw` paths).
#'
#' @param sample_ids A `character` vector specifying the sample IDs
#' to include in the plot. If `NULL`, the function uses the sample IDs
#' specified in the `lcmsPlot` object.
#' @param type A `character` value indicating the type of plot;
#' one of `"boxplot"`, `"violin"`, `"jitter"`.
#' @return This function returns another function that takes an `lcmsPlot`
#' object and produces a modified version containing the generated
#' total ion current (TIC) in its `data` slot. It is designed to be used
#' with the `+` operator, which serves as a layering mechanism.
#' Each use of `+` incrementally enriches the `lcmsPlot` object by
#' adding new data or visual components.
#' @export
#' @examples
#' data_obj <- get_XCMSnExp_object_example()
#'
#' p <- lcmsPlot(data_obj, sample_id_column = "sample_name") +
#'   lp_total_ion_current(type = "violin") +
#'   lp_arrange(group_by = "sample_id")
#' p
lp_total_ion_current <- function(sample_ids = NULL, type = "boxplot") {
    function(obj) {
        if (!is_xcms_data(obj@data@data_obj) && !is.character(obj@data@data_obj)) {
            stop("total_ion_current: to plot the total ion current the data object should be of class XCMSnExp, MsExperiment, or a character vector of raw file paths.")
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
#' The `lp_intensity_map` function produces an intensity map
#' in which signal intensity is represented at each m/z / RT coordinate.
#'
#' @param mz_range A `numeric` value indicating the m/z range of the map.
#' @param rt_range A `numeric` value indicating the RT range of the map.
#' @param sample_ids A `character` vector specifying the sample IDs
#' to include in the plot. If `NULL`, the function uses the sample IDs
#' specified in the `lcmsPlot` object.
#' @param x_dim A `character` value indicating which dimension to place on the
#' x-axis. One of `"rt"` (default) or `"mz"`.
#' @param y_dim A `character` value indicating which dimension to place on the
#' y-axis. One of `"mz"` (default) or `"rt"`.
#' @param fill_scale A ggplot2 scale object to use for the fill aesthetic
#' (e.g. `scale_fill_gradient(low = "white", high = "red")`). When `NULL`
#' (default), uses `scale_fill_viridis_c` for intensity and
#' `scale_fill_viridis_d` for density plots. Applies to `geom = "tile"` and
#' `geom = "density"`.
#' @param geom A `character` value selecting how the map is drawn. One of
#' `"tile"` (default) for a binned heatmap, `"point"` for a scatter of the
#' individual centroids, or `"density"` for a
#' 2D kernel density estimate. `"point"` skips the binning entirely, so gaps in
#' the mass traces stay visible rather than being implied away.
#' @param point_size A `numeric` value controlling the point size when
#' `geom = "point"`.
#' @param bin_rt A `numeric` value giving the retention-time bin width in
#' seconds used when `geom = "tile"`. Defaults to `0.1`, which is too fine for
#' slow-scanning data and too coarse for fast-scanning data.
#' @param bin_mz A `numeric` value giving the m/z bin width used when
#' `geom = "tile"`. Defaults to `0.1`; high-resolution data usually wants less.
#' @param colour_scale A ggplot2 scale object to use for the colour aesthetic
#' when `geom = "point"`. When `NULL` (default), uses `scale_colour_viridis_c`.
#' @return This function returns another function that takes an
#' `lcmsPlot` object and produces a modified version containing the generated
#' 2D intensity map in its `data` slot. It is designed to be used with the
#' `+` operator, which serves as a layering mechanism.
#' Each use of `+` incrementally enriches the `lcmsPlot` object by
#' adding new data or visual components.
#' @export
#' @examples
#' raw_files <- dir(
#'   system.file("cdf", package = "faahKO"),
#'   full.names = TRUE,
#'   recursive = TRUE)[1]
#'
#' p <- lcmsPlot(raw_files) +
#'   lp_intensity_map(
#'     mz_range = c(200, 600),
#'     rt_range = c(4200, 4500),
#'     geom = "density")
#' p
lp_intensity_map <- function(
    mz_range,
    rt_range,
    sample_ids = NULL,
    x_dim = "rt",
    y_dim = "mz",
    fill_scale = NULL,
    geom = "tile",
    point_size = 0.5,
    bin_rt = 0.1,
    bin_mz = 0.1,
    colour_scale = NULL
) {
    geom <- match.arg(geom, c("tile", "point", "density"))
    x_dim <- match.arg(x_dim, c("rt", "mz"))
    y_dim <- match.arg(y_dim, c("mz", "rt"))

    function(obj) {
        if (is.null(sample_ids)) {
            sample_ids <- obj@data@metadata$sample_id
        }

        obj@options$intensity_maps <- list(
            show = TRUE,
            sample_ids = sample_ids,
            mz_range = mz_range,
            rt_range = rt_range,
            x_dim = x_dim,
            y_dim = y_dim,
            fill_scale = fill_scale,
            geom = geom,
            point_size = point_size,
            bin_rt = bin_rt,
            bin_mz = bin_mz,
            colour_scale = colour_scale
        )

        obj@data <- create_intensity_map(obj@data, obj@options)

        return(obj)
    }
}

#' Generate the retention time difference plot between raw and adjusted datasets
#'
#' The `lp_rt_diff_plot` function generates the data necessary to plot
#' the difference between the raw and retention time adjusted datasets.
#' Only applicable to `XCMSnExp` and `MsExperiment` objects.
#'
#' @return This function returns another function that takes an
#' `lcmsPlot` object and produces a modified version containing the generated
#' retention time differences, between raw and adjusted,  in its `data` slot.
#' It is designed to be used with the `+` operator, which serves as a layering
#' mechanism. Each use of `+` incrementally enriches the `lcmsPlot` object by
#' adding new data or visual components.
#' @export
#' @examples
#' data_obj <- get_XCMSnExp_object_example(
#'   indices = 1:3,
#'   should_group_peaks = TRUE)
#' p <- lcmsPlot(data_obj, sample_id_column = "sample_name") +
#'   lp_rt_diff_plot()
#' p
lp_rt_diff_plot <- function() {
    function(obj) {
        if (!is_xcms_data(obj@data@data_obj)) {
            stop("lp_rt_diff_plot: to plot the RT differences the data object should be either of class XCMSnExp or MsExperiment.")
        }

        if (!xcms_utils$has_rt_alignment_been_performed(obj@data@data_obj)) {
            stop("lp_rt_diff_plot: RT alignment was not performed.")
        }

        obj@options$rt_diff <- list(show = TRUE)
        obj@data <- create_rt_diff(obj@data, obj@options)

        return(obj)
    }
}

#' Define the arrangement of chromatograms
#'
#' The `lp_arrange` function specifies how chromatograms should be arranged
#' when visualised. It determines the grouping metadata factor through
#' the `group_by` parameter.
#'
#' @param group_by A `character` value determining the
#' column to group by in the samples metadata.
#' @return A function that takes an `lcmsPlot` object and returns a modified
#' version with the specified arrangement options
#' stored in `options$arrangement`. It is intended for use with the `+`
#' operator, which incrementally layers new data or visual components
#' onto the `lcmsPlot` object.
#' @export
#' @examples
#' raw_files <- dir(
#'    system.file("cdf", package = "faahKO"),
#'    full.names = TRUE,
#'    recursive = TRUE)[1:5]
#'
#' ## Plots chromatograms overlayed without specifying a grouping factor
#' p <- lcmsPlot(raw_files) +
#'   lp_chromatogram(aggregation_fun = "max")
#' p
#'
#' ## Plots chromatograms overlayed specifying a grouping factor
#' ## (e.g., sample_id)
#' p <- p + lp_arrange(group_by = "sample_id")
#' p
lp_arrange <- function(group_by) {
    make_interface_function(
        name = "lp_arrange",
        args_list = as.list(environment()),
        fn = function(obj) {
            obj@options$arrangement <- list(
                group_by = group_by
            )
            return(obj)
        }
    )
}

#' Define the plot's faceting
#'
#' The `lp_facets` function arranges plots into a grid based
#' on a metadata factor, creating a series of smaller plots (facets).
#'
#' @param facets A `character` vector of factors from the sample metadata to use
#' for faceting.
#' @param ncol A `numeric` value indicating the number of columns in the layout.
#' @param nrow A `numeric` value indicating the number of rows in the layout.
#' @param free_x A `logical` value indicating whether the x-axis scales
#' are allowed to vary across panels.
#' @param free_y A `logical` value indicating whether the y-axis scales
#' are allowed to vary across panels.
#' @return A function that takes an `lcmsPlot` object and returns a modified
#' version with the specified faceting options stored in `options$facets`.
#' It is intended for use with the `+` operator, which incrementally layers
#' new data or visual components onto the `lcmsPlot` object.
#' @export
#' @examples
#' raw_files <- dir(
#'    system.file("cdf", package = "faahKO"),
#'    full.names = TRUE,
#'    recursive = TRUE)[1:5]
#'
#' ## Plots chromatograms overlayed
#' p <- lcmsPlot(raw_files) +
#'   lp_chromatogram(aggregation_fun = "max")
#' p
#'
#' ## Using lp_facets we create facets for each sample_id
#' p <- p + lp_facets(facets = "sample_id")
#' p
lp_facets <- function(
    facets,
    ncol = NULL,
    nrow = NULL,
    free_x = FALSE,
    free_y = FALSE
) {
    make_interface_function(
        name = "lp_facets",
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
#' The `lp_grid` function arranges plots into a matrix of panels defined
#' by row and column faceting metadata factors.
#'
#' @param rows A `character` value indicating the factors that
#' represent rows.
#' @param cols A `character` value indicating the factors that
#' represent columns.
#' @param free_x A `logical` value indicating whether the x-axis scales
#' are allowed to vary across panels.
#' @param free_y A `logical` value indicating whether the y-axis scales
#' are allowed to vary across panels.
#' @return A function that takes an `lcmsPlot` object and returns a modified
#' version with the specified grid options stored in `options$grid`.
#' It is intended for use with the `+` operator, which incrementally layers
#' new data or visual components onto the `lcmsPlot` object.
#' @export
#' @examples
#' raw_files <- dir(
#'   system.file("cdf", package = "faahKO"),
#'   full.names = TRUE,
#'   recursive = TRUE
#' )[1:4]
#'
#' ## Create metadata for the samples
#' metadata <- data.frame(
#'   sample_id = sub("\\.CDF", "", basename(raw_files)),
#'   factor1 = c("S", "S", "C", "C"),
#'   factor2 = c("T", "U", "T", "U")
#' )
#'
#' ## Create feature chromatograms for the specified samples
#' p <- lcmsPlot(raw_files, metadata = metadata) +
#'   lp_chromatogram(features = rbind(c(
#'     mzmin = 334.9,
#'     mzmax = 335.1,
#'     rtmin = 2700,
#'     rtmax = 2900)))
#' p
#'
#' ## Arrange chromatograms in a grid split by experimental factors
#' ## Rows correspond to `factor1` and columns correspond to `factor2`
#' p <- p + lp_grid(rows = "factor1", cols = "factor2")
#' p
lp_grid <- function(rows, cols, free_x = FALSE, free_y = FALSE) {
    make_interface_function(
        name = "lp_grid",
        args_list = as.list(environment()),
        fn = function(obj) {
            obj@options$grid <- list(
                rows = rows,
                cols = cols,
                free_x = free_x,
                free_y = free_y
            )
            return(obj)
        }
    )
}

#' Define the labels of the plot
#'
#' The `lp_labels` function allows the specification
#' of the plot title and the legend title.
#'
#' @param title A `character` value indicating the plot title.
#' @param legend A `character` value indicating the legend's title.
#' @return A function that takes an `lcmsPlot` object and returns a modified
#' version with the specified label options stored in `options$labels`.
#' It is intended for use with the `+` operator, which incrementally layers
#' new data or visual components onto the `lcmsPlot` object.
#' @export
#' @examples
#' raw_files <- dir(
#'    system.file("cdf", package = "faahKO"),
#'    full.names = TRUE,
#'    recursive = TRUE)[1:5]
#'
#' ## Create a chromatogram plot by grouping samples into batches
#' ## By default, the legend is derived from the grouping variable
#' p <- lcmsPlot(raw_files, batch_size = 2) +
#'   lp_chromatogram(features = rbind(c(
#'     mzmin = 334.9,
#'     mzmax = 335.1,
#'     rtmin = 2700,
#'     rtmax = 2900))) +
#'   lp_arrange(group_by = "sample_id")
#' p
#'
#' ## Customise the legend label
#' p <- p + lp_labels(legend = "Sample")
#' p
lp_labels <- function(title = NULL, legend = NULL) {
    make_interface_function(
        name = "lp_labels",
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
#' @param position A `character` value indicating the legend's position.
#' One of `"top"`, `"right"`, `"bottom"`, `"left"`, or `"inside"`.
#' @return A function that takes an `lcmsPlot` object and returns a modified
#' version with the specified legend options stored in `options$legend`.
#' It is intended for use with the `+` operator, which incrementally layers
#' new data or visual components onto the `lcmsPlot` object.
#' @export
#' @examples
#' raw_files <- dir(
#'    system.file("cdf", package = "faahKO"),
#'    full.names = TRUE,
#'    recursive = TRUE)[1:5]
#'
#' ## Create a chromatogram plot grouped by sample with a custom legend label
#' p <- lcmsPlot(raw_files, batch_size = 2) +
#'   lp_chromatogram(features = rbind(c(
#'     mzmin = 334.9,
#'     mzmax = 335.1,
#'     rtmin = 2700,
#'     rtmax = 2900))) +
#'   lp_arrange(group_by = "sample_id") +
#'   lp_labels(legend = "Sample")
#' p
#'
#' ## Move the legend below the plot
#' p <- p + lp_legend(position = "bottom")
#' p
lp_legend <- function(position = NULL)  {
    make_interface_function(
        name = "lp_legend",
        args_list = as.list(environment()),
        fn = function(obj) {
            obj@options$legend <- list(
                position = position
            )
            return(obj)
        }
    )
}

#' Define a vertical line on a retention time value
#'
#' @param intercept A `numeric` value indicating the retention time
#' axis (x-axis) intercept.
#' @param line_type A `character` value indicating the line type.
#' One of `"solid"`, `"dashed"`, `"dotted"`,
#' `"dotdash"`, `"longdash"`, `"twodash"`.
#' @param color A `character` value indicating the line color.
#' @return A function that takes an `lcmsPlot` object and returns a modified
#' version with the specified RT line options stored in `options$rt_lines`.
#' It is intended for use with the `+` operator, which incrementally layers
#' new data or visual components onto the `lcmsPlot` object.
#' @export
#' @examples
#' raw_files <- dir(
#'    system.file("cdf", package = "faahKO"),
#'    full.names = TRUE,
#'    recursive = TRUE)[1:4]
#'
#' ## Create chromatogram plots faceted by sample
#' p <- lcmsPlot(raw_files) +
#'   lp_chromatogram(features = rbind(c(
#'     mzmin = 334.9,
#'     mzmax = 335.1,
#'     rtmin = 2700,
#'     rtmax = 2900))) +
#'   lp_facets(facets = 'sample_id', ncol = 4)
#' p
#'
#' ## Add a vertical retention time reference line
#' p <- p + lp_rt_line(intercept = 2800, line_type = 'solid', color = 'red')
#' p
lp_rt_line <- function(intercept, line_type = 'dashed', color = 'black') {
    make_interface_function(
        name = "lp_rt_line",
        args_list = as.list(environment()),
        fn = function(obj) {
            rt_line_obj <- list(
                intercept = intercept,
                line_type = line_type,
                color = color
            )
            obj@options$rt_lines <- append(
                obj@options$rt_lines,
                list(rt_line_obj))
            return(obj)
        }
    )
}

#' Define the plot layout
#'
#' @param design Specification of the location of areas in the layout
#' See https://patchwork.data-imaginist.com/reference/wrap_plots.html
#' @return A function that takes an `lcmsPlot` object and returns a modified
#' version with the specified layout options stored in `options$layout`.
#' It is intended for use with the `+` operator, which incrementally layers
#' new data or visual components onto the `lcmsPlot` object.
#' @export
#' @examples
#' data_obj <- get_XCMSnExp_object_example(indices = 1)
#'
#' ## Plot chromatograms and spectra for selected samples and features
#' p <- lcmsPlot(data_obj, sample_id_column = 'sample_name') +
#'   lp_chromatogram(
#'     features = rbind(
#'       c(mzmin = 334.9, mzmax = 335.1, rtmin = 2700, rtmax = 2900),
#'       c(mzmin = 278.99721, mzmax = 279.00279, rtmin = 2740, rtmax = 2840)
#'     ),
#'     sample_ids = 'ko15',
#'     highlight_peaks = TRUE
#'   ) +
#'   lp_spectra(mode = "closest_apex", ms_level = 1) +
#'   lp_facets(facets = "feature_id", ncol = 2)
#'
#' ## Customise panel layout to place chromatogram above spectra
#' p <- p + lp_layout(design = "C\nS\nS")
lp_layout <- function(design = NULL) {
    make_interface_function(
        name = "lp_layout",
        args_list = as.list(environment()),
        fn = function(obj) {
            obj@options$layout <- list(
                design = design
            )
            return(obj)
        }
    )
}

#' Define the options to use when plotting LC-MS data coming from
#' Compound Discoverer results.
#'
#' @param compounds_query A `character` value indicating the expression
#' used to filter compounds from the Compound Discoverer results.
#' The expression is evaluated on the compound table and can reference
#' the following columns:
#' \describe{
#'   \item{name}{Compound name.}
#'   \item{formula}{Chemical formula of the compound.}
#'   \item{adduct}{Ion adduct (e.g. `[M+H]+`, `[M-H]-`).}
#'   \item{rt}{Retention time of the compound (in seconds).}
#'   \item{rtmin}{Minimum retention time of the compound peak.}
#'   \item{rtmax}{Maximum retention time of the compound peak.}
#'   \item{mz}{Mass-to-charge ratio (m/z) of the detected ion.}
#'   \item{maxo}{Maximum observed peak intensity.}
#'   \item{into}{Integrated peak area reported by Compound Discoverer.}
#' }
#' @param rt_extend A `numeric` value indicating how much (in seconds)
#' the retention time window should be extended on each side of the
#' compound peak when extracting and plotting chromatograms.
#' @return A function that takes an `lcmsPlot` object and returns a modified
#' version with the specified Compound Discoverer options stored in
#' `options$compound_discoverer`. It is intended for use with the `+` operator,
#' which incrementally layers new data or visual components onto
#' the `lcmsPlot` object.
#' @export
#' @examples
#' ## A minimal stand-in for a Compound Discoverer Scripting Node export.
#' ## See [CompoundDiscovererNodeSource()] for the table layout.
#' node_dir <- tempfile("cd_node")
#' dir.create(node_dir)
#'
#' write_tab <- function(df, name) {
#'     path <- file.path(node_dir, name)
#'     utils::write.table(
#'         df, path, sep = "\t", row.names = FALSE, quote = FALSE)
#'     path
#' }
#'
#' node_args <- list(Tables = list(
#'     list(DataFile = write_tab(
#'         data.frame(
#'             "Compounds ID" = c(1, 2),
#'             "Name" = c("L-Proline", "L-Kynurenine"),
#'             "Formula" = c("C5 H9 N O2", "C10 H12 N2 O3"),
#'             check.names = FALSE),
#'         "compounds.txt")),
#'     list(DataFile = write_tab(
#'         data.frame(
#'             "Compounds per File ID" = c(1, 2),
#'             "StudyFileID" = c(1, 2),
#'             "Area" = c(24050, 18700),
#'             "Intensity" = c(24050, 18700),
#'             check.names = FALSE),
#'         "cpf.txt")),
#'     list(DataFile = write_tab(
#'         data.frame(
#'             "Features ID" = c(1, 2),
#'             "mz" = c(116.0704, 209.0918),
#'             "Ion" = c("[M+H]+1", "[M+H]+1"),
#'             "Area" = c(24050, 18700),
#'             "RT [min]" = c(7.10, 6.10),
#'             "Left RT [min]" = c(6.95, 5.95),
#'             "Right RT [min]" = c(7.25, 6.25),
#'             check.names = FALSE),
#'         "features.txt")),
#'     list(DataFile = write_tab(
#'         data.frame(
#'             "Compounds ID" = c(1, 2),
#'             "Compounds per File ID" = c(1, 2),
#'             check.names = FALSE),
#'         "link_cmp_cpf.txt")),
#'     list(DataFile = write_tab(
#'         data.frame(
#'             "Compounds per File ID" = c(1, 2),
#'             "Features ID" = c(1, 2),
#'             check.names = FALSE),
#'         "link_cpf_feat.txt")),
#'     list(DataFile = write_tab(
#'         data.frame(
#'             "StudyFileID" = c(1, 2),
#'             "File Name" = c("l-proline-MS1.mzML", "l-kynurenine-MS1.mzML"),
#'             check.names = FALSE),
#'         "study_files.txt"))))
#'
#' json_path <- file.path(node_dir, "node_args.json")
#' writeLines(jsonlite::toJSON(node_args, auto_unbox = TRUE), json_path)
#'
#' mzml_dir <- file.path(tempdir(), "cd-node-example")
#' utils::unzip(
#'     system.file("extdata", "standards-mzml.zip", package = "lcmsPlot"),
#'     exdir = mzml_dir)
#'
#' ds <- CompoundDiscovererNodeSource(
#'     node_args = json_path,
#'     sample_paths = file.path(
#'         mzml_dir, "mzml",
#'         c("l-proline-MS1.mzML", "l-kynurenine-MS1.mzML")))
#'
#' ## `compounds_query` is evaluated on the compound table.
#' p <- lcmsPlot(ds) +
#'   lp_compound_discoverer(
#'     compounds_query = 'name %in% c("L-Proline", "L-Kynurenine")',
#'     rt_extend = 15
#'   ) +
#'   lp_chromatogram(highlight_peaks = TRUE) +
#'   lp_grid(rows = "sample_id", cols = "name", free_x = TRUE) +
#'   lp_labels(title = "Compound Discoverer example", legend = "Sample") +
#'   lp_legend(position = "bottom")
#' p
lp_compound_discoverer <- function(compounds_query = NULL, rt_extend = 10) {
    make_interface_function(
        name = "lp_compound_discoverer",
        args_list = as.list(environment()),
        fn = function(obj) {
            if (!is_cd_result(obj@data@data_obj) &&
                !is_cd_node_source(obj@data@data_obj)) {
                stop("lp_compound_discoverer: The data object is not a Compound Discoverer results connection or scripting-node source")
            }

            obj@options$compound_discoverer <- list(
                compounds_query = compounds_query,
                rt_extend = rt_extend
            )

            return(obj)
        }
    )
}

#' Define the options to use when plotting LC-MS data coming from
#' a LipidSearch result file.
#'
#' @param lipids_query A `character` value indicating the expression used to
#' select which lipids to plot. The expression is evaluated on the lipid table
#' and can reference the following columns:
#' \describe{
#'   \item{name}{Plot label: the lipid ion, disambiguated with its retention
#'   time when the same ion is reported more than once.}
#'   \item{lipid_ion}{Lipid ion as reported by LipidSearch, e.g. `AEA(18:2)+H`
#'   (4.2 `LipidIon` / 5.2 `LipidID`).}
#'   \item{lipid_group}{LipidSearch lipid group.}
#'   \item{class}{Lipid class, e.g. `PC`, `AcCa`.}
#'   \item{sub_class}{Lipid subclass (5.2 only; `NA` for 4.2).}
#'   \item{fatty_acid}{Fatty acid composition (4.2 only; `NA` for 5.2).}
#'   \item{formula}{Ion formula.}
#'   \item{adduct}{Adduct: the explicit `AdductIon` column for 5.2 (e.g. `M+H`),
#'   otherwise derived from the ion name for 4.2 (e.g. `+H`).}
#'   \item{calc_mz}{Theoretical m/z.}
#'   \item{rej}{LipidSearch rejection flag as a `logical` (`TRUE` if rejected).}
#'   \item{mz}{m/z used for extraction (observed where available).}
#'   \item{rt, rtmin, rtmax}{Retention time and peak bounds, in seconds.}
#'   \item{into}{Integrated peak area; `NA` where nothing was detected.}
#'   \item{maxo}{Peak height; `NA` where nothing was detected.}
#'   \item{grade}{LipidSearch identification grade (`A` to `D`).}
#'   \item{mscore, sn}{LipidSearch m-score / ID score and signal-to-noise
#'   (`sn` is 4.2 only).}
#'   \item{detected}{Whether LipidSearch integrated a peak in this sample.}
#'   \item{lipid_rank}{Dense rank of lipids by descending total area, for
#'   top-N selection.}
#' }
#'
#' The query selects *lipids*, not lipid/sample rows: a predicate such as
#' `grade == "A"` keeps the lipids graded A in at least one sample, and each of
#' them is then extracted from **every** sample - including samples that the
#' result file does not cover, which use the lipid's consensus m/z and
#' retention-time window.
#' @param rt_extend A `numeric` value indicating how much (in seconds) the
#' retention time window should be extended on each side of the lipid peak when
#' extracting and plotting chromatograms.
#' @return A function that takes an `lcmsPlot` object and returns a modified
#' version with the specified LipidSearch options stored in
#' `options$lipid_search`. It is intended for use with the `+` operator,
#' which incrementally layers new data or visual components onto
#' the `lcmsPlot` object.
#' @seealso [LipidSearchSource()]
#' @export
#' @examples
#' ## A minimal LipidSearch 4.2 export; see [LipidSearchSource()] for the
#' ## format details.
#' results_path <- file.path(tempdir(), "lipids_pos.txt")
#'
#' columns <- c(
#'     "Rej.", "LipidIon", "LipidGroup", "Class", "FattyAcid", "CalcMz",
#'     "IonFormula",
#'     "Area[c-1]", "Area[c-2]", "Height[c-1]", "Height[c-2]",
#'     "Rt[c-1]", "Rt[c-2]", "ObsMz[c-1]", "ObsMz[c-2]",
#'     "Grade[c-1]", "Grade[c-2]")
#'
#' rows <- list(
#'     c(0, "PRO(0:0)+H", "PRO(0:0)", "PRO", "(0:0)", 116.0706,
#'       "C5 H10 O2 N1", 24050, 0, 24050, 0, 7.10, 0,
#'       "116.0704", "", "A", ""),
#'     c(0, "KYN(0:0)+H", "KYN(0:0)", "KYN", "(0:0)", 209.0921,
#'       "C10 H13 O3 N2", 0, 18700, 0, 18700, 0, 6.10,
#'       "", "209.0918", "", "A"))
#'
#' writeLines(
#'     c("#[c-1]:l-proline-MS1.raw",
#'       "#[c-2]:l-kynurenine-MS1.raw",
#'       "#mScoreThreshold:5.0",
#'       "",
#'       vapply(
#'           c(list(columns), rows),
#'           function(x) paste0(paste(x, collapse = "\t"), "\t"),
#'           character(1))),
#'     results_path)
#'
#' mzml_dir <- file.path(tempdir(), "lipid-search-example")
#' utils::unzip(
#'     system.file("extdata", "standards-mzml.zip", package = "lcmsPlot"),
#'     exdir = mzml_dir)
#'
#' ds <- LipidSearchSource(
#'     results_path = results_path,
#'     sample_paths = file.path(
#'         mzml_dir, "mzml",
#'         c("l-proline-MS1.mzML", "l-kynurenine-MS1.mzML")))
#'
#' ## `lipids_query` selects lipids, which are then extracted from every sample.
#' p <- lcmsPlot(ds) +
#'   lp_lipid_search(
#'     lipids_query = 'grade == "A"',
#'     rt_extend = 30
#'   ) +
#'   lp_chromatogram(highlight_peaks = TRUE) +
#'   lp_grid(rows = "sample_id", cols = "name", free_x = TRUE) +
#'   lp_labels(title = "LipidSearch example", legend = "Sample") +
#'   lp_legend(position = "bottom")
#' p
lp_lipid_search <- function(lipids_query = NULL, rt_extend = 30) {
    make_interface_function(
        name = "lp_lipid_search",
        args_list = as.list(environment()),
        fn = function(obj) {
            if (!is_lipid_search_source(obj@data@data_obj)) {
                stop("lp_lipid_search: The data object is not a LipidSearch source")
            }

            obj@options$lipid_search <- list(
                lipids_query = lipids_query,
                rt_extend = rt_extend
            )

            return(obj)
        }
    )
}

#' Plot chromatographic peak boundaries in the rt / m/z plane
#'
#' `lp_chrom_peak_rects()` draws one rectangle per detected chromatographic
#' peak, spanning `rtmin`-`rtmax` by `mzmin`-`mzmax`. It works in two ways,
#' matching the two xcms methods that draw the same primitive.
#'
#' Called on its own it owns a panel, drawing the rectangles on an otherwise
#' empty retention time / m/z frame as `xcms::plotChromPeaks()` does, with one
#' facet per sample. Called after `lp_mass_trace()` or `lp_intensity_map()` it
#' instead decorates that panel, as `xcms::plot(type = "XIC")` does, and follows
#' the host's axes, window and samples. On an intensity-versus-rt panel the
#' equivalent is
#' `lp_chromatogram(highlight_peaks = TRUE, highlight_peaks_mode = "rectangle")`.
#'
#' The data object has to report chromatographic peaks. On m/z-binned data
#' `mzmin` equals `mzmax`, so each box is given a minimum height and reads as a
#' horizontal segment rather than vanishing.
#'
#' @param sample_ids A `character` vector of sample IDs to include. `NULL`
#' (default) follows the host layer when there is one, so the overlay never
#' introduces panels the host does not draw, and otherwise uses every sample.
#' @param border A `character` colour for the rectangle outline.
#' @param fill A `character` colour for the rectangle interior, or `NA`
#' (default) for unfilled boxes.
#' @param alpha A `numeric` value in `[0, 1]` controlling opacity.
#' @param linewidth A `numeric` value controlling the outline width.
#' @param rt_range An optional length-2 `numeric` restricting the retention-time
#' range, equivalent to `xlim` in `xcms::plotChromPeaks()`. Ignored when the
#' layer decorates an `lp_intensity_map()` panel, which supplies its own window.
#' @param mz_range An optional length-2 `numeric` restricting the m/z range,
#' equivalent to `ylim` in `xcms::plotChromPeaks()`. Ignored when the layer
#' decorates an `lp_intensity_map()` panel.
#' @return A layer function for use with the `+` operator.
#' @export
#' @examples
#' data_obj <- get_XCMSnExp_object_example(
#'   indices = 1:2,
#'   should_group_peaks = TRUE)
#'
#' ## On its own the layer owns the panel, one facet per sample.
#' p <- lcmsPlot(data_obj, sample_id_column = "sample_name") +
#'   lp_chrom_peak_rects(rt_range = c(2500, 3500), mz_range = c(300, 320))
#' p
#'
#' ## After a host layer it decorates that panel instead.
#' p <- lcmsPlot(data_obj, sample_id_column = "sample_name") +
#'   lp_intensity_map(mz_range = c(300, 320), rt_range = c(2500, 3500)) +
#'   lp_chrom_peak_rects()
#' p
lp_chrom_peak_rects <- function(
        sample_ids = NULL,
        border = "#c0392b",
        fill = NA,
        alpha = 0.25,
        linewidth = 0.3,
        rt_range = NULL,
        mz_range = NULL
) {
    make_interface_function(
        name = "lp_chrom_peak_rects",
        args_list = as.list(environment()),
        fn = function(obj) {
            has_map <- isTRUE(obj@options$intensity_maps$show)
            has_trace <- isTRUE(obj@options$mass_traces$show)

            if (has_trace && !has_map && nrow(obj@data@mass_traces) == 0) {
                stop(
                    "lp_chrom_peak_rects: the mass trace panel is empty, so ",
                    "there is nothing to overlay. Mass traces are built from ",
                    "the raw files, which XChromatograms and XChromatogram ",
                    "objects do not carry; use lp_intensity_map() as the host ",
                    "panel instead."
                )
            }

            # lp_intensity_map() does not extract peaks, so unlike the
            # chromatogram layers this overlay may have to populate them.
            if (nrow(obj@data@detected_peaks) == 0) {
                peaks <- get_detected_peaks(obj@data@data_obj)

                if (is.null(peaks) || nrow(peaks) == 0) {
                    stop(
                        "lp_chrom_peak_rects: the data object has no detected ",
                        "peaks. Provide an XCMSnExp or MsExperiment object ",
                        "with chromatographic peaks."
                    )
                }

                obj@data@detected_peaks <- peaks |>
                    left_join(obj@data@metadata, by = "sample_index")
                validObject(obj@data)
            }

            # Fail here rather than drawing nothing: an empty host panel is
            # dropped before the renderer runs, so a missing boundary column
            # would otherwise surface as a silently absent overlay.
            missing_cols <- setdiff(
                c("rtmin", "rtmax", "mzmin", "mzmax"),
                names(obj@data@detected_peaks))
            if (length(missing_cols) > 0) {
                stop(
                    "lp_chrom_peak_rects: the detected peaks are missing the ",
                    "column(s) ", paste(missing_cols, collapse = ", "),
                    ", so the peak boundaries cannot be drawn."
                )
            }

            # Follow the host panel's samples by default, or the overlay would
            # contribute peaks for samples the host does not draw and faceting
            # would invent empty panels for them. Standalone there is no host to
            # follow, so every sample is drawn.
            if (is.null(sample_ids)) {
                sample_ids <- if (has_map) {
                    obj@options$intensity_maps$sample_ids
                } else if (has_trace) {
                    obj@options$chromatograms$sample_ids
                } else {
                    NULL
                }
            }

            obj@options$chrom_peak_rects <- list(
                show = TRUE,
                sample_ids = sample_ids,
                border = border,
                fill = fill,
                alpha = alpha,
                linewidth = linewidth,
                rt_range = rt_range,
                mz_range = mz_range
            )

            return(obj)
        }
    )
}

#' Overlay precursor ion purity scores on a chromatogram
#'
#' `lp_purity_overlay()` adds a `geom_point` layer to the chromatogram panel
#' where each triangle marks the retention time of an MS/MS fragmentation event,
#' coloured by the interpolated precursor ion purity (`inPurity`). Requires the
#' data object to be a `purityA` result from the msPurity package and
#' `lp_chromatogram()` to be called first.
#'
#' @param sample_ids A `character` vector of sample IDs to include.
#' `NULL` uses all samples.
#' @param threshold A `numeric` value in `[0, 1]` drawn as a reference
#' on the colour scale midpoint. `NULL` defaults to `0.5`.
#' @param point_size A `numeric` value controlling the point size.
#' @return A layer function for use with the `+` operator.
#' @export
#' @examplesIf requireNamespace("msPurity", quietly = TRUE)
#' ## lcmsPlot ships two DDA files carrying real precursor isolation metadata.
#' mzml_dir <- file.path(tempdir(), "purity-example")
#' utils::unzip(
#'     system.file("extdata", "standards-mzml.zip", package = "lcmsPlot"),
#'     exdir = mzml_dir)
#'
#' ms2_files <- file.path(
#'     mzml_dir, "mzml",
#'     c("l-proline-MS2.mzML", "l-kynurenine-MS2.mzML"))
#'
#' pa <- msPurity::purityA(ms2_files)
#'
#' ## The overlay annotates an existing chromatogram, so the retention time
#' ## window has to span the fragmentation events to show anything.
#' features <- data.frame(
#'     sample_id = c("l-proline-MS2", "l-kynurenine-MS2"),
#'     mz = c(116.0703793, 209.091824),
#'     rt = c(425.74, 365.72))
#'
#' p <- lcmsPlot(pa) +
#'     lp_chromatogram(features = features, ppm = 10, rt_tol = 60) +
#'     lp_purity_overlay(threshold = 0.7)
#' p
lp_purity_overlay <- function(
        sample_ids = NULL,
        threshold = NULL,
        point_size = 2
) {
    make_interface_function(
        name = "lp_purity_overlay",
        args_list = as.list(environment()),
        fn = function(obj) {
            if (!is_mspurity_data(obj@data@data_obj)) {
                stop("lp_purity_overlay requires a purityA data object")
            }

            if (!isTRUE(obj@options$chromatograms$show)) {
                stop("lp_purity_overlay must be called after lp_chromatogram()")
            }

            obj@options$purity_overlay <- list(
                show = TRUE,
                sample_ids = sample_ids,
                threshold = threshold,
                point_size = point_size
            )
            obj@options$purity_scores_sample_ids <- sample_ids

            obj@data <- create_purity_scores(obj@data, obj@options)
            return(obj)
        }
    )
}

#' Visualise the isolation window for a precursor fragmentation event
#'
#' `lp_isolation_window()` extracts the MS1 survey spectrum at the scan
#' immediately preceding a chosen MS2 event from a `purityA` object and
#' annotates it with a semi-transparent rectangle spanning the isolation window
#' and a dashed line at the precursor m/z. The `inPurity` score is shown as a
#' plot subtitle.
#'
#' @param pid An integer or integer vector of precursor IDs (the `pid` column in
#' `purityA@@puritydf`) selecting which event(s) to display. `NULL` shows all
#' events (use with care on large datasets).
#' @param sample_id A `character` sample ID to filter by when `pid` is `NULL`.
#' @param half_width A `numeric` value (in Da) for the isolation window
#' half-width on each side of the precursor m/z. Default is `0.5`.
#' @param zoom_factor A `numeric` multiplier controlling how far beyond the
#' isolation window the x-axis extends. The visible range is
#' `precursor_mz ± zoom_factor * half_width`. Default is `3`.
#' @return A layer function for use with the `+` operator.
#' @export
#' @examplesIf requireNamespace("msPurity", quietly = TRUE)
#' mzml_dir <- file.path(tempdir(), "purity-example")
#' utils::unzip(
#'     system.file("extdata", "standards-mzml.zip", package = "lcmsPlot"),
#'     exdir = mzml_dir)
#'
#' ms2_files <- file.path(
#'     mzml_dir, "mzml",
#'     c("l-proline-MS2.mzML", "l-kynurenine-MS2.mzML"))
#'
#' pa <- msPurity::purityA(ms2_files)
#'
#' ## Inspect the fragmentation events available to plot.
#' slot(pa, "puritydf")[, c("pid", "precursorMZ", "precursorRT", "inPurity")]
#'
#' ## `pid = 4` is a co-isolated precursor (inPurity ~ 0.51). These files
#' ## record an isolation width of 0.6 Da on each side of the precursor.
#' p <- lcmsPlot(pa) +
#'     lp_isolation_window(pid = 4, half_width = 0.6)
#' p
lp_isolation_window <- function(
        pid = NULL,
        sample_id = NULL,
        half_width = 0.5,
        zoom_factor = 3
) {
    make_interface_function(
        name = "lp_isolation_window",
        args_list = as.list(environment()),
        fn = function(obj) {
            if (!is_mspurity_data(obj@data@data_obj)) {
                stop("lp_isolation_window requires a purityA data object")
            }

            obj@options$isolation_window <- list(
                show = TRUE,
                half_width = half_width,
                zoom_factor = zoom_factor
            )
            obj@options$spectra$show <- TRUE

            obj@data <- create_isolation_window(
                obj@data, obj@options, pid, sample_id)
            return(obj)
        }
    )
}

#' Plot precursor ion purity scores as a timeline
#'
#' `lp_purity_timeline()` generates a scatter plot of `inPurity` (y-axis)
#' versus retention time (x-axis), one point per MS/MS acquisition event,
#' coloured by sample. An optional horizontal dashed line marks a purity
#' threshold. Requires the data object to be a `purityA` result.
#'
#' @param sample_ids A `character` vector of sample IDs to include.
#' `NULL` uses all samples.
#' @param threshold A `numeric` value in `[0, 1]` drawn as a horizontal
#' reference line. `NULL` suppresses the line.
#' @return A layer function for use with the `+` operator.
#' @export
#' @examplesIf requireNamespace("msPurity", quietly = TRUE)
#' mzml_dir <- file.path(tempdir(), "purity-example")
#' utils::unzip(
#'     system.file("extdata", "standards-mzml.zip", package = "lcmsPlot"),
#'     exdir = mzml_dir)
#'
#' ms2_files <- file.path(
#'     mzml_dir, "mzml",
#'     c("l-proline-MS2.mzML", "l-kynurenine-MS2.mzML"))
#'
#' pa <- msPurity::purityA(ms2_files)
#'
#' ## One point per MS/MS event, coloured by sample.
#' p <- lcmsPlot(pa) +
#'     lp_purity_timeline(threshold = 0.7)
#' p
lp_purity_timeline <- function(
        sample_ids = NULL,
        threshold = NULL
) {
    make_interface_function(
        name = "lp_purity_timeline",
        args_list = as.list(environment()),
        fn = function(obj) {
            if (!is_mspurity_data(obj@data@data_obj)) {
                stop("lp_purity_timeline requires a purityA data object")
            }

            obj@options$purity_timeline <- list(
                show = TRUE,
                sample_ids = sample_ids,
                threshold = threshold
            )
            obj@options$purity_scores_sample_ids <- sample_ids

            obj@data <- create_purity_scores(obj@data, obj@options)
            return(obj)
        }
    )
}

#' Plot the distribution of precursor ion purity scores per sample
#'
#' `lp_purity_distribution()` generates a violin (or box/jitter) plot of
#' `inPurity` scores grouped by sample. An optional horizontal dashed line
#' marks a purity acceptance threshold. Requires the data object to be a
#' `purityA` result from the msPurity package.
#'
#' @param sample_ids A `character` vector of sample IDs to include.
#' `NULL` uses all samples.
#' @param threshold A `numeric` value in `[0, 1]` drawn as a horizontal
#' reference line. `NULL` suppresses the line.
#' @param type A `character` value; one of `"violin"`, `"boxplot"`,
#' or `"jitter"`. Defaults to `"violin"`.
#' @return A layer function for use with the `+` operator.
#' @export
#' @examplesIf requireNamespace("msPurity", quietly = TRUE)
#' mzml_dir <- file.path(tempdir(), "purity-example")
#' utils::unzip(
#'     system.file("extdata", "standards-mzml.zip", package = "lcmsPlot"),
#'     exdir = mzml_dir)
#'
#' ms2_files <- file.path(
#'     mzml_dir, "mzml",
#'     c("l-proline-MS2.mzML", "l-kynurenine-MS2.mzML"))
#'
#' pa <- msPurity::purityA(ms2_files)
#'
#' ## Compare the spread of purity scores between samples.
#' p <- lcmsPlot(pa) +
#'     lp_purity_distribution(threshold = 0.7, type = "boxplot")
#' p
lp_purity_distribution <- function(
        sample_ids = NULL,
        threshold = NULL,
        type = "violin"
) {
    make_interface_function(
        name = "lp_purity_distribution",
        args_list = as.list(environment()),
        fn = function(obj) {
            if (!is_mspurity_data(obj@data@data_obj)) {
                stop("lp_purity_distribution requires a purityA data object")
            }

            obj@options$purity_distribution <- list(
                show = TRUE,
                sample_ids = sample_ids,
                threshold = threshold,
                type = type
            )
            obj@options$purity_scores_sample_ids <- sample_ids

            obj@data <- create_purity_scores(obj@data, obj@options)
            return(obj)
        }
    )
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
