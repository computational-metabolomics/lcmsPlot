.validators <- list(
    metadata = function(df) {
        validate_data_frame(df, list(
            field("sample_index", is.numeric),
            field("sample_id", is.character),
            field("sample_path", is.character, required = FALSE)
        ))
    },
    feature_metadata = function(df) {
        validate_data_frame(df, list(
            field("feature_metadata_id", is.numeric)
        ))
    },
    chromatograms = function(df) {
        validate_data_frame(df, list(
            field("rt", is.numeric),
            field("intensity", is.numeric),
            field("metadata_index", is.numeric),
            field("feature_metadata_id", is.numeric)
        ))
    },
    mass_traces = function(df) {
        validate_data_frame(df, list(
            field("rt", is.numeric),
            field("mz", is.numeric),
            field("metadata_index", is.numeric),
            field("feature_metadata_id", is.numeric)
        ), exact = TRUE)
    },
    spectra = function(df) {
        validate_data_frame(df, list(
            field("mz", is.numeric),
            field("intensity", is.numeric),
            field("rt", is.numeric),
            field("metadata_index", is.numeric),
            field("feature_metadata_id", is.numeric),
            field("reference", is.logical)
        ), exact = TRUE)
    },
    total_ion_current = function(df) {
        validate_data_frame(df, list(
            field("intensity", is.numeric),
            field("metadata_index", is.numeric),
            field("feature_metadata_id", is.numeric)
        ), exact = TRUE)
    },
    intensity_maps = function(df) {
        validate_data_frame(df, list(
            field("rt", is.numeric),
            field("mz", is.numeric),
            field("intensity", is.numeric),
            field("metadata_index", is.numeric),
            field("feature_metadata_id", is.numeric)
        ), exact = TRUE)
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
            field("metadata_index", is.numeric),
            field("feature_metadata_id", is.numeric)
        ))
    },
    rt_diff = function(df) {
        validate_data_frame(df, list(
            field("rt_raw", is.numeric),
            field("rt_adj", is.numeric),
            field("diff", is.numeric),
            field("metadata_index", is.numeric),
            field("feature_metadata_id", is.numeric)
        ), exact = TRUE)
    },
    detected_peaks = function(df) {
        validate_data_frame(df, list(
            field("mz", is.numeric),
            field("rt", is.numeric),
            field("rtmin", is.numeric),
            field("rtmax", is.numeric),
            field("sample_index", is.numeric),
            field("sample_id", is.character),
            field("sample_path", is.character, required = FALSE)
        ))
    },
    purity_scores = function(df) {
        validate_data_frame(df, list(
            field("rt", is.numeric),
            field("in_purity", is.numeric),
            field("precursor_mz", is.numeric),
            field("metadata_index", is.numeric),
            field("feature_metadata_id", is.numeric)
        ), exact = TRUE)
    }
)

DATASET_TYPES <- c(
    "chromatograms",
    "mass_traces",
    "spectra",
    "peak_density",
    "intensity_maps",
    "total_ion_current",
    "rt_diff"
)

#' Create an instance of class `lcmsPlotDataContainer` from a data object
#'
#' The `create_data_container_from_obj` function creates an instance
#' of class `lcmsPlotDataContainer` given a data object.
#' See `lcmsPlotDataContainer` for more information about the
#' supported data objects.
#'
#' @param data_obj The data object (see `lcmsPlotDataContainer`).
#' @param sample_id_column A `character` value indicating the sample ID column.
#' @param metadata A `data.frame` containing the samples metadata
#' in case it is not provided in the dataset object.
#' @return An instance of class `lcmsPlotDataContainer`. The object
#' contains the input data and the standardised metadata.
#' @export
#' @examples
#' raw_files <- dir(
#'    system.file("cdf", package = "faahKO"),
#'    full.names = TRUE,
#'    recursive = TRUE)[1:5]
#'
#' data_container <- create_data_container_from_obj(
#'   data_obj = raw_files,
#'   sample_id_column = NULL,
#'   metadata = NULL
#' )
#' @keywords internal
create_data_container_from_obj <- function(
    data_obj,
    sample_id_column,
    metadata
) {
    if (is_cd_results_path(data_obj)) {
        data_obj <- open_cd_result_connection(data_obj)
    }

    new("lcmsPlotDataContainer",
        data_obj = data_obj,
        metadata = get_metadata(data_obj, sample_id_column, metadata),
        chromatograms = tibble(),
        mass_traces = tibble(),
        spectra = tibble(),
        peak_density = tibble(),
        total_ion_current = tibble(),
        intensity_maps = tibble(),
        rt_diff = tibble(),
        feature_metadata = tibble(
            feature_metadata_id = numeric(),
            metadata_index = numeric()
        ),
        detected_peaks = tibble(),
        purity_scores = tibble())
}

#' A unified storing mechanism for LC-MS data
#'
#' The `lcmsPlotDataContainer` class allows the storage of different
#' types of LC-MS data.
#' This class can be used independently from the plotting utilities,
#' however the preferred approach is to use it with the `lcmsPlotClass` class.
#'
#' @slot data_obj The data object. One of: `XCMSnExp`, `MsExperiment`,
#' `MChromatograms`, `XChromatograms`, `XChromatogram`, `XcmsRawList`,
#' `purityA`, or `character` representing mzML paths.
#' @slot metadata A `data.frame` containing the sample metadata.
#' @slot chromatograms A `data.frame` containing the chromatograms.
#' @slot mass_traces A `data.frame` containing the mass traces.
#' @slot spectra A `data.frame` containing the spectra.
#' @slot peak_density A `data.frame` containing peak density curve data and
#' optional feature-group rectangles, as produced by `lp_peak_density()`.
#' @slot total_ion_current A `data.frame` containing the total ion current.
#' @slot intensity_maps A `data.frame` containing the 2D intensity maps
#' representing the distribution of detected peaks across m/z and RT.
#' @slot rt_diff A `data.frame` containing the raw and adjusted RT values.
#' @slot feature_metadata A `data.frame` containing feature/compound annotations
#' attached to datasets through a column called `feature_metadata_id`.
#' @slot detected_peaks A `data.frame` containing the detected peaks from
#' an `XCMSnExp`, `MsExperiment`, or `purityA` object.
#' @slot purity_scores A `data.frame` containing per-scan precursor ion purity
#' scores from a `purityA` object, populated by the msPurity layer functions.
#' @export
setClass(
    "lcmsPlotDataContainer",
    slots = list(
        data_obj = "ANY",
        metadata = "data.frame",
        chromatograms = "data.frame",
        mass_traces = "data.frame",
        spectra = "data.frame",
        peak_density = "data.frame",
        total_ion_current = "data.frame",
        intensity_maps = "data.frame",
        rt_diff = "data.frame",
        feature_metadata = "data.frame",
        detected_peaks = "data.frame",
        purity_scores = "data.frame"
    ),
    prototype = list(
        data_obj = NULL,
        metadata = NULL,
        chromatograms = NULL,
        mass_traces = NULL,
        spectra = NULL,
        peak_density = NULL,
        total_ion_current = NULL,
        intensity_maps = NULL,
        rt_diff = NULL,
        feature_metadata = NULL,
        detected_peaks = NULL,
        purity_scores = NULL
    )
)

setValidity("lcmsPlotDataContainer", function(object) {
    ret <- TRUE

    obj_types <- c(
        "XCMSnExp",
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
        "character")

    if (!inherits(object@data_obj, obj_types)) {
        ret <- sprintf(
            "@data_obj must inherit from one of: %s",
            paste(obj_types, collapse = ", ")
        )
    } else {
        ret <- validate_object(object, .validators)
    }

    ret
})

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
#' data_obj <- new("lcmsPlotDataContainer",
#'     data_obj = raw_files,
#'     metadata = tibble::tibble(),
#'     chromatograms = tibble::tibble(),
#'     mass_traces = tibble::tibble(),
#'     spectra = tibble::tibble(),
#'     peak_density = tibble::tibble(),
#'     total_ion_current = tibble::tibble(),
#'     intensity_maps = tibble::tibble(),
#'     rt_diff = tibble::tibble(),
#'     feature_metadata = tibble::tibble(),
#'     detected_peaks = tibble::tibble(),
#'     purity_scores = tibble::tibble())
#' data_obj
setMethod(
    f = "show",
    signature = "lcmsPlotDataContainer",
    function(object) {
        cat("Object of class", class(object), "\n")
        cat(" Data object type:", class(object@data_obj), "\n")

        print_df_dim <- function(x, name) {
            if (is.null(x)) {
                cat(" ", name, ": NULL\n")
            } else if (is.data.frame(x)) {
                d <- dim(x)
                cat(" ", name, ":", d[1], "rows x", d[2], "columns\n")
            } else {
                cat(" ", name, ": not a data frame\n")
            }
        }

        print_df_dim(object@metadata, "metadata")
        print_df_dim(object@chromatograms, "chromatograms")
        print_df_dim(object@mass_traces, "mass_traces")
        print_df_dim(object@spectra, "spectra")
        print_df_dim(object@peak_density, "peak_density")
        print_df_dim(object@total_ion_current, "total_ion_current")
        print_df_dim(object@intensity_maps, "intensity_maps")
        print_df_dim(object@rt_diff, "rt_diff")
        print_df_dim(object@feature_metadata, "feature_metadata")
        print_df_dim(object@detected_peaks, "detected_peaks")
        print_df_dim(object@purity_scores, "purity_scores")
    }
)
