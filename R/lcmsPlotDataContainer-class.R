#' A unified storing mechanism for LC-MS data
#'
#' The `lcmsPlotDataContainer` class holds the input data and every dataset
#' derived from it for plotting, using the RforMassSpectrometry data structures
#' where one exists.
#'
#' Every dataset is keyed by `sample_index`, the row of the sample metadata
#' (`sampleMetadata()`) it belongs to. Feature annotations (such as `feature_id`, `name`, or `formula`)
#' are carried by the rows they describe: `chromData()` columns of the
#' chromatograms, spectra variables of the spectra, or columns of the
#' summaries. Sample metadata is joined in once, when the plot is rendered.
#'
#' @slot experiment An `MsExperiment` whose `sampleData()` holds the sample
#' metadata, with columns `sample_index` (the row number), `sample_id`, and
#' `sample_path`. Its spectra are only loaded once a layer needs the raw data.
#' @slot source The input object, kept for the information that has no
#' RforMassSpectrometry counterpart (detected peaks, feature definitions,
#' compound tables, or purity scores). One of: `XcmsExperiment`,
#' `MsExperiment`, `MChromatograms`, `XChromatograms`, `XChromatogram`,
#' `XcmsRawList`, `ExternalDataSource`, `CompoundDiscovererNodeSource`,
#' `LipidSearchSource`, `DBIConnection`, `purityA`, or `character` (raw file
#' paths). `XCMSnExp` input is stored as an `XcmsExperiment`.
#' @slot chromatograms A `Chromatograms` object with the extracted
#' chromatograms, held in memory.
#' @slot spectra A `Spectra` object with the extracted spectra, held in memory.
#' Spectra variables `sample_index`, `feature_id`, and `reference` link each
#' spectrum to its sample, feature, and library match.
#' @slot detected_peaks A `data.frame` of the detected peaks selected by the
#' chromatogram layer, with columns `mz`, `mzmin`, `mzmax`, `rt`, `rtmin`,
#' `rtmax`, `into`, `maxo`, and `sample_index`, plus `feature_id` when the peak
#' was matched to a feature.
#' @slot summaries A named `list` of `data.frame`s holding the derived datasets
#' that have no RforMassSpectrometry counterpart: `mass_traces`,
#' `intensity_maps`, `total_ion_current`, `rt_diff`, `peak_density`,
#' `peak_count_image`, and `purity_scores`.
#' @export
setClass(
    "lcmsPlotDataContainer",
    slots = list(
        experiment = "MsExperiment",
        source = "ANY",
        chromatograms = "Chromatograms",
        spectra = "Spectra",
        detected_peaks = "data.frame",
        summaries = "list"
    ),
    prototype = list(
        experiment = MsExperiment::MsExperiment(),
        source = NULL,
        chromatograms = Chromatograms::Chromatograms(),
        spectra = Spectra::Spectra(),
        detected_peaks = tibble(),
        summaries = list()
    )
)

.source_types <- c(
    "XcmsExperiment",
    "MsExperiment",
    "MChromatograms",
    "XChromatograms",
    "XChromatogram",
    "XcmsRawList",
    "ExternalDataSource",
    "CompoundDiscovererNodeSource",
    "LipidSearchSource",
    "DBIConnection",
    "purityA",
    "character"
)

.summary_validators <- list(
    mass_traces = function(df) {
        validate_data_frame(df, list(
            field("rt", is.numeric),
            field("mz", is.numeric),
            field("sample_index", is.numeric)
        ))
    },
    intensity_maps = function(df) {
        validate_data_frame(df, list(
            field("rt", is.numeric),
            field("mz", is.numeric),
            field("intensity", is.numeric),
            field("sample_index", is.numeric)
        ))
    },
    total_ion_current = function(df) {
        validate_data_frame(df, list(
            field("intensity", is.numeric),
            field("sample_index", is.numeric)
        ))
    },
    rt_diff = function(df) {
        validate_data_frame(df, list(
            field("rt_raw", is.numeric),
            field("rt_adj", is.numeric),
            field("diff", is.numeric),
            field("sample_index", is.numeric)
        ))
    },
    peak_density = function(df) {
        validate_data_frame(df, list(
            field("rt", is.numeric),
            field("density", is.numeric),
            field("rtmin", is.numeric),
            field("rtmax", is.numeric),
            field("data_type", is.character),
            field("mzmin", is.numeric),
            field("mzmax", is.numeric),
            field("feature_index", is.numeric),
            field("sample_index", is.numeric)
        ))
    },
    peak_count_image = function(df) {
        validate_data_frame(df, list(
            field("rt", is.numeric),
            field("n_peaks", is.numeric),
            field("sample_index", is.numeric)
        ))
    },
    purity_scores = function(df) {
        validate_data_frame(df, list(
            field("rt", is.numeric),
            field("in_purity", is.numeric),
            field("precursor_mz", is.numeric),
            field("sample_index", is.numeric)
        ))
    }
)

setValidity("lcmsPlotDataContainer", function(object) {
    if (!is.null(object@source) && !inherits(object@source, .source_types)) {
        return(sprintf(
            "@source must inherit from one of: %s",
            paste(.source_types, collapse = ", ")))
    }

    sample_data <- MsExperiment::sampleData(object@experiment)
    if (nrow(sample_data) > 0) {
        missing_cols <- setdiff(
            c("sample_index", "sample_id"), colnames(sample_data))
        if (length(missing_cols) > 0) {
            return(paste(
                "the sample metadata is missing the column(s):",
                paste(missing_cols, collapse = ", ")))
        }
        if (!identical(
            as.integer(sample_data$sample_index), seq_len(nrow(sample_data)))) {
            return("the sample metadata's sample_index must be the row number")
        }
    }

    unknown <- setdiff(names(object@summaries), names(.summary_validators))
    if (length(unknown) > 0) {
        return(paste(
            "Unknown summaries:", paste(unknown, collapse = ", ")))
    }

    for (name in names(object@summaries)) {
        result <- .summary_validators[[name]](object@summaries[[name]])
        if (!isTRUE(result)) {
            return(paste0("summaries$", name, ": ", result))
        }
    }

    TRUE
})

#' Create an instance of class `lcmsPlotDataContainer` from a data object
#'
#' The `lcmsPlotDataContainer` function standardises the input data object
#' into an `lcmsPlotDataContainer`. See `lcmsPlotDataContainer-class` for the
#' supported data objects.
#'
#' Creating the container is cheap: raw spectra are only read once a layer
#' needs them.
#'
#' @param data_obj The data object (see `lcmsPlotDataContainer-class`).
#' @param sample_id_column A `character` value indicating the sample ID column.
#' @param metadata A `data.frame` containing the samples metadata
#' in case it is not provided in the dataset object.
#' @return An instance of class `lcmsPlotDataContainer`.
#' @export
#' @examples
#' raw_files <- dir(
#'    system.file("cdf", package = "faahKO"),
#'    full.names = TRUE,
#'    recursive = TRUE)[1:5]
#'
#' data_container <- lcmsPlotDataContainer(
#'   data_obj = raw_files,
#'   sample_id_column = NULL,
#'   metadata = NULL
#' )
#' sampleMetadata(data_container)
lcmsPlotDataContainer <- function(
    data_obj,
    sample_id_column = "sample_id",
    metadata = NULL
) {
    if (is_cd_results_path(data_obj)) {
        data_obj <- open_cd_result_connection(data_obj)
    }

    if (inherits(data_obj, "XCMSnExp")) {
        data_obj <- methods::as(data_obj, "XcmsExperiment")
    }

    sample_metadata <- get_metadata(data_obj, sample_id_column, metadata) |>
        mutate(sample_index = row_number())

    if (inherits(data_obj, "MsExperiment")) {
        experiment <- data_obj
        MsExperiment::sampleData(experiment) <-
            S4Vectors::DataFrame(sample_metadata, check.names = FALSE)
    } else {
        experiment <- MsExperiment::MsExperiment(
            sampleData = S4Vectors::DataFrame(
                sample_metadata, check.names = FALSE))
    }

    if (is(data_obj, "XcmsRawList")) {
        sps <- xcmsraw_list_to_spectra(data_obj)
        experiment <- .link_spectra(experiment, sps, sps$sample_index)
    }

    obj <- new(
        "lcmsPlotDataContainer",
        experiment = experiment,
        source = data_obj)
    validObject(obj)
    obj
}

#' Show a summary of an instance of class `lcmsPlotDataContainer`
#'
#' @param object An instance of class `lcmsPlotDataContainer`.
#' @return Invisible \code{NULL}
#' @export
#' @examples
#' raw_files <- dir(
#'    system.file("cdf", package = "faahKO"),
#'    full.names = TRUE,
#'    recursive = TRUE)[1:5]
#'
#' lcmsPlotDataContainer(raw_files)
setMethod(
    f = "show",
    signature = "lcmsPlotDataContainer",
    function(object) {
        sample_data <- MsExperiment::sampleData(object@experiment)
        cat("Object of class", class(object), "\n")
        cat(" Source type:", class(object@source)[1], "\n")
        cat(" Samples:", nrow(sample_data), "\n")
        cat(" Raw spectra loaded:",
            length(MsExperiment::spectra(object@experiment)) > 0, "\n")
        cat(" Chromatograms:", length(object@chromatograms), "\n")
        cat(" Spectra:", length(object@spectra), "\n")
        cat(" Detected peaks:", nrow(object@detected_peaks), "\n")
        for (name in names(object@summaries)) {
            cat(" ", name, ":", nrow(object@summaries[[name]]), "rows\n")
        }
        invisible(NULL)
    }
)

#' Access the data held by an `lcmsPlotDataContainer`
#'
#' Accessors for the data structures held by an `lcmsPlotDataContainer`, or by
#' the container of an `lcmsPlotClass` object.
#'
#' - `sampleMetadata()` returns the sample metadata as a `tibble`, with
#'   `sample_index`, `sample_id`, and `sample_path` columns. It is the
#'   `sampleData()` of the underlying `MsExperiment`.
#' - `chromatograms()` returns the extracted chromatograms as a
#'   `Chromatograms` object.
#' - `spectra()` returns the extracted spectra as a `Spectra` object.
#' - `detectedPeaks()` returns the detected peaks selected by the chromatogram
#'   layer as a `data.frame`.
#' - `summaries()` returns the derived datasets as a named `list` of
#'   `data.frame`s.
#'
#' @param object An `lcmsPlotDataContainer` or `lcmsPlotClass` object.
#' @param ... Ignored.
#' @return See the description.
#' @name lcmsPlotDataContainer-accessors
#' @examples
#' raw_files <- dir(
#'    system.file("cdf", package = "faahKO"),
#'    full.names = TRUE,
#'    recursive = TRUE)[1:2]
#'
#' p <- lcmsPlot(raw_files) + lp_chromatogram()
#' sampleMetadata(p)
#' chromatograms(p)
NULL

#' @rdname lcmsPlotDataContainer-accessors
#' @export
setGeneric("sampleMetadata", function(object, ...) {
    standardGeneric("sampleMetadata")
})

#' @rdname lcmsPlotDataContainer-accessors
#' @export
setGeneric("detectedPeaks", function(object, ...) {
    standardGeneric("detectedPeaks")
})

#' @rdname lcmsPlotDataContainer-accessors
#' @export
setGeneric("summaries", function(object, ...) {
    standardGeneric("summaries")
})

#' @rdname lcmsPlotDataContainer-accessors
#' @export
setMethod("sampleMetadata", "lcmsPlotDataContainer", function(object, ...) {
    .sample_table(object)
})

#' @rdname lcmsPlotDataContainer-accessors
#' @importFrom ProtGenerics spectra
#' @export
setMethod("spectra", "lcmsPlotDataContainer", function(object, ...) {
    object@spectra
})

#' @rdname lcmsPlotDataContainer-accessors
#' @importFrom ProtGenerics chromatograms
#' @export
setMethod("chromatograms", "lcmsPlotDataContainer", function(object, ...) {
    object@chromatograms
})

#' @rdname lcmsPlotDataContainer-accessors
#' @export
setMethod("detectedPeaks", "lcmsPlotDataContainer", function(object, ...) {
    object@detected_peaks
})

#' @rdname lcmsPlotDataContainer-accessors
#' @export
setMethod("summaries", "lcmsPlotDataContainer", function(object, ...) {
    object@summaries
})

# The sample metadata as a tibble, the form every join uses.
.sample_table <- function(container) {
    as_tibble(as.data.frame(
        MsExperiment::sampleData(container@experiment),
        optional = TRUE))
}

# Join the sample metadata onto a dataset keyed by `sample_index`. Columns the
# dataset already carries win, so annotations are never suffixed .x/.y.
.join_sample_data <- function(df, sample_table) {
    if (!"sample_index" %in% colnames(df)) {
        return(df)
    }
    clashing <- setdiff(intersect(colnames(df), colnames(sample_table)),
                        "sample_index")
    sample_table <- sample_table[, setdiff(colnames(sample_table), clashing),
                                 drop = FALSE]
    df$sample_index <- as.integer(df$sample_index)
    left_join(as_tibble(df), sample_table, by = "sample_index")
}

# Restrict the sample table to the requested sample IDs, keeping its order.
.select_samples <- function(container, sample_ids = NULL) {
    samples <- .sample_table(container)
    if (is.null(sample_ids)) {
        return(samples)
    }
    samples |> filter(.data$sample_id %in% sample_ids)
}

.set_summary <- function(container, name, df) {
    container@summaries[[name]] <- df
    validObject(container)
    container
}

.summary <- function(container, name) {
    df <- container@summaries[[name]]
    if (is.null(df)) tibble() else df
}
