.validators <- list(
    additional_metadata = function(df) {
        TRUE
    },
    chromatograms = function(df) {
        cols <- c(
            "rt", "intensity", "metadata_index", "additional_metadata_index")
        nrow(df) == 0 || identical(colnames(df), cols)
    },
    mass_traces = function(df) {
        cols <- c("rt", "mz", "metadata_index", "additional_metadata_index")
        nrow(df) == 0 || identical(colnames(df), cols)
    },
    spectra = function(df) {
        cols <- c(
            "mz", "intensity", "rt", "metadata_index",
            "additional_metadata_index", "reference")
        nrow(df) == 0 || identical(colnames(df), cols)
    },
    total_ion_current = function(df) {
        cols <- c("intensity", "metadata_index", "additional_metadata_index")
        nrow(df) == 0 || identical(colnames(df), cols)
    },
    intensity_maps = function(df) {
        cols <- c(
            "rt", "mz", "intensity",
            "metadata_index", "additional_metadata_index")
        nrow(df) == 0 || identical(colnames(df), cols)
    },
    rt_diff = function(df) {
        cols <- c(
            "rt_raw", "rt_adj", "diff",
            "metadata_index", "additional_metadata_index")
        nrow(df) == 0 || identical(colnames(df), cols)
    },
    detected_peaks = function(df) {
        cols <- c("mz", "rt", "rtmin", "rtmax", "sample_index")
        nrow(df) == 0 || all(cols %in% colnames(df))
    }
)

DATASET_TYPES <- c(
    "chromatograms",
    "mass_traces",
    "spectra",
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
        chromatograms = data.frame(),
        mass_traces = data.frame(),
        spectra = data.frame(),
        total_ion_current = data.frame(),
        intensity_maps = data.frame(),
        rt_diff = data.frame(),
        additional_metadata = data.frame(),
        detected_peaks = data.frame())
}

#' A unified storing mechanism for LC-MS data
#'
#' The `lcmsPlotDataContainer` class allows the storage of different
#' types of LC-MS data.
#' This class can be used independently from the plotting utilities,
#' however the preferred approach is to use it with the `lcmsPlotClass` class.
#'
#' @slot data_obj The data object. One of: `XCMSnExp`, `MsExperiment` or
#' `character` representing mzML paths.
#' @slot metadata A `data.frame` containing the sample metadata.
#' @slot chromatograms A `data.frame` containing the chromatograms.
#' @slot mass_traces A `data.frame` containing the mass traces.
#' @slot spectra A `data.frame` containing the spectra.
#' @slot total_ion_current A `data.frame` containing the total ion current.
#' @slot intensity_maps A `data.frame` containing the 2D intensity maps
#' representing the distribution of detected peaks across m/z and RT.
#' @slot rt_diff A `data.frame` containing the raw and adjusted RT values.
#' @slot additional_metadata A `data.frame` containing additional information
#' attached to datasets through a column called `additional_metadata_index`.
#' @slot detected_peaks A `data.frame` containing the detected peaks from
#' an `XCMSnExp` or `MsExperiment` object.
#' @export
setClass(
    "lcmsPlotDataContainer",
    slots = list(
        data_obj = "ANY",
        metadata = "data.frame",
        chromatograms = "data.frame",
        mass_traces = "data.frame",
        spectra = "data.frame",
        total_ion_current = "data.frame",
        intensity_maps = "data.frame",
        rt_diff = "data.frame",
        additional_metadata = "data.frame",
        detected_peaks = "data.frame"
    ),
    prototype = list(
        data_obj = NULL,
        metadata = NULL,
        chromatograms = NULL,
        mass_traces = NULL,
        spectra = NULL,
        total_ion_current = NULL,
        intensity_maps = NULL,
        rt_diff = NULL,
        additional_metadata = NULL,
        detected_peaks = NULL
    )
)

setValidity("lcmsPlotDataContainer", function(object) {
    ret <- TRUE

    obj_types <- c(
        "XCMSnExp",
        "MsExperiment",
        "ExternalDataSource",
        "DBIConnection",
        "character")

    if (!inherits(object@data_obj, obj_types)) {
        ret <- sprintf(
            "@data_obj must inherit from one of: %s",
            paste(obj_types, collapse = ", ")
        )
    } else {
        for (validator_name in names(.validators)) {
            df <- slot(object, validator_name)

            if (!.validators[[validator_name]](df)) {
                ret <- paste0(validator_name, " did not pass validation.")
                break
            }
        }
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
#'     metadata = data.frame(),
#'     chromatograms = data.frame(),
#'     mass_traces = data.frame(),
#'     spectra = data.frame(),
#'     total_ion_current = data.frame(),
#'     intensity_maps = data.frame(),
#'     rt_diff = data.frame(),
#'     additional_metadata = data.frame(),
#'     detected_peaks = data.frame())
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
                cat(" ", name, ": not a data.frame\n")
            }
        }

        print_df_dim(object@metadata, "metadata")
        print_df_dim(object@chromatograms, "chromatograms")
        print_df_dim(object@mass_traces, "mass_traces")
        print_df_dim(object@spectra, "spectra")
        print_df_dim(object@total_ion_current, "total_ion_current")
        print_df_dim(object@intensity_maps, "intensity_maps")
        print_df_dim(object@rt_diff, "rt_diff")
        print_df_dim(object@additional_metadata, "additional_metadata")
        print_df_dim(object@detected_peaks, "detected_peaks")
    }
)
