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
    BPPARAM = BiocParallel::SerialParam()
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

            if (nrow(data_df) == 0) {
                stop("Empty dataset ", dataset_name)
            }

            data_df <- merge_by_index(
                data_df,
                object@data@metadata,
                index_col = 'metadata_index'
            )

            if (nrow(object@data@additional_metadata) > 0) {
                data_df <- merge_by_index(
                    data_df,
                    object@data@additional_metadata,
                    index_col = 'additional_metadata_index'
                )
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

.has_data <- function(object) {
    data_slots <- c(
        "chromatograms",
        "mass_traces",
        "spectra",
        "total_ion_current",
        "intensity_maps",
        "rt_diff",
        "additional_metadata",
        "detected_peaks"
    )

    any(vapply(
        data_slots,
        function(s) nrow(slot(object@data, s)),
        integer(1)
    ) > 0)
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
                object <- .render_plot(object, additional_datasets = list())
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
#' in an LC–MS dataset are either summed to produce the total ion current (TIC)
#' chromatogram or the most intense peak is selected to produce the
#' base peak chromatogram (BPC).
#' To create such chromatograms do not specify the `features` parameter as
#' that will create the chromatograms for the selected features.
#' In this context, the main parameter if `aggregation_fun` which can take
#' either `sum` (TIC) or `max` (BPC).
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
#' @param highlight_peaks A `logical` value indicating whether to highlight
#' the detected peaks; the input data must be an `XCMSnExp` or `MsExperiment`
#' object.
#' @param highlight_peaks_color A `character` value indicating the color of the
#' highlighted peaks.
#' @param highlight_peaks_factor A `character` value indicating the factor from
#' the metadata that determines the color. By default it colors by `sample_id`.
#' @param aggregation_fun A `character` value indicating which aggregation
#' function to use for the spectra intensities; one of `max` or `sum`.
#' Only applicable to summary chromatograms.
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
#' @param highlight_apices A `logical` value indicating whether to
#' highlight apices with the corresponding RT values in a chromatogram.
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
    highlight_peaks = FALSE,
    highlight_peaks_color = NULL,
    highlight_peaks_factor = "sample_id",
    aggregation_fun = "max",
    rt_type = "uncorrected",
    rt_unit = "second",
    intensity_unit = "absolute",
    fill_gaps = FALSE,
    highlight_apices = list(column = NULL, top_n = NULL)
) {
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
                rt_tol = rt_tol,
                highlight_peaks = highlight_peaks,
                highlight_peaks_color = highlight_peaks_color,
                highlight_peaks_factor = highlight_peaks_factor,
                aggregation_fun = aggregation_fun,
                rt_type = rt_type,
                rt_unit = rt_unit,
                intensity_unit = intensity_unit,
                fill_gaps = fill_gaps,
                highlight_apices = highlight_apices
            )

            if (is.null(features)) {
                if (is_cd_result(obj@data@data_obj)) {
                    obj@data <- create_chromatograms_from_compound_discoverer(
                        obj@data,
                        obj@options)
                } else {
                    obj@data <- create_full_rt_chromatograms(
                        obj@data,
                        obj@options)
                }
            } else if (is.character(features)) {
                obj@data <- create_chromatograms_from_feature_ids(
                    obj@data,
                    obj@options)
            } else {
                obj@data <- create_chromatograms_from_features(
                    obj@data,
                    obj@options)
            }

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
    match_target_index = NULL
) {
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
                match_target_index = match_target_index
            )

            obj@data <- create_spectra(obj@data, obj@options)
            return(obj)
        }
    )
}

#' Define the total ion current (TIC)
#'
#' The `lp_total_ion_current` generates summary data for the
#' total ion current (TIC) of the selected samples.
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
        if (!is_xcms_data(obj@data@data_obj)) {
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
#' The `lp_intensity_map` function produces an intensity map
#' in which retention time is shown on the x-axis, m/z on the y-axis,
#' and signal intensity is represented at each corresponding coordinate.
#'
#' @param mz_range A `numeric` value indicating the m/z range of the map.
#' @param rt_range A `numeric` value indicating the RT range of the map.
#' @param sample_ids A `character` vector specifying the sample IDs
#' to include in the plot. If `NULL`, the function uses the sample IDs
#' specified in the `lcmsPlot` object.
#' @param density A `logical` value indicating whether to show a density plot.
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
#'     density = TRUE)
#' p
lp_intensity_map <- function(
    mz_range,
    rt_range,
    sample_ids = NULL,
    density = FALSE
) {
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
#' \dontrun{
#' lcmsPlot("cd_example.cdResult") +
#'   lp_compound_discoverer(
#'     compounds_query = 'name %in% c("Proline", "Betaine")',
#'     rt_extend = 5
#'   ) +
#'   lp_chromatogram(highlight_peaks = TRUE) +
#'   lp_grid(rows = "sample_id", cols = "name", free_x = TRUE) +
#'   lp_labels(title = "Compound Discoverer example", legend = "Sample") +
#'   lp_legend(position = "bottom")
#' }
lp_compound_discoverer <- function(compounds_query = NULL, rt_extend = 10) {
    make_interface_function(
        name = "lp_compound_discoverer",
        args_list = as.list(environment()),
        fn = function(obj) {
            if (!is_cd_result(obj@data@data_obj)) {
                stop("lp_compound_discoverer: The data object is not a Compound Discoverer DB connection")
            }

            obj@options$compound_discoverer <- list(
                compounds_query = compounds_query,
                rt_extend = rt_extend
            )

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
        obj <- .render_plot(obj, additional_datasets = list())
        return(obj@plot)
    }
}
