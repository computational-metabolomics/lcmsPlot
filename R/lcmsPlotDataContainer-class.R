.validators <- list(
  additional_metadata = function(df) {
    nrow(df) == 0 || "metadata_index" %in% colnames(df)
  },
  chromatograms = function(df) {
    nrow(df) == 0 || identical(colnames(df), c("rt", "intensity", "metadata_index", "additional_metadata_index"))
  },
  mass_traces = function(df) {
    nrow(df) == 0 || identical(colnames(df), c("rt", "mz", "metadata_index", "additional_metadata_index"))
  },
  spectra = function(df) {
    nrow(df) == 0 || identical(colnames(df), c("mz", "intensity", "rt", "metadata_index", "additional_metadata_index", "reference"))
  },
  total_ion_current = function(df) {
    nrow(df) == 0 || identical(colnames(df), c("intensity", "metadata_index", "additional_metadata_index"))
  },
  intensity_maps = function(df) {
    nrow(df) == 0 || identical(colnames(df), c("rt", "mz", "intensity", "metadata_index", "additional_metadata_index"))
  },
  rt_diff = function(df) {
    nrow(df) == 0 || identical(colnames(df), c("rt_raw", "rt_adj", "diff", "metadata_index", "additional_metadata_index"))
  },
  detected_peaks = function(df) {
    nrow(df) == 0 || all(c("mz", "rt", "rtmin", "rtmax", "sample_index") %in% colnames(df))
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

#' Create an lcmsPlotDataContainer object from a data object (e.g., XCMSnExp).
#'
#' @param data_obj The data object (e.g., XCMSnExp).
#' @param sample_id_column The sample ID column.
#' @param metadata The sample metadata in case it's not provided in the data object.
#' @returns The created lcmsPlotDataContainer object.
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
create_data_container_from_obj <- function(data_obj, sample_id_column, metadata) {
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

#' lcmsPlotDataContainer class.
#'
#' @slot data_obj The data object. One of:
#' - XCMSnExp
#' - MsExperiment
#' - character: vector of sample paths
#' @slot metadata The sample metadata.
#' @slot chromatograms The chromatograms.
#' @slot mass_traces The mass traces.
#' @slot spectra The spectra.
#' @slot total_ion_current The Total Ion Current (TIC).
#' @slot intensity_maps The 2D intensity maps.
#' @slot rt_diff The RT difference dataset for RT alignment.
#' @slot additional_metadata Additional information attached to datasets.
#' @slot detected_peaks The detected peaks.
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

  if (!inherits(object@data_obj, c("XCMSnExp", "MsExperiment", "character"))) {
    ret <- "@data_obj must be either XCMSnExp, MsExperiment, or character."
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

#' Check if the lcmsPlotDataContainer data object is an xcms object.
#'
#' @param object The lcmsPlotDataContainer object.
#' @returns Whether the data object is an xcms object.
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
#'
#' is_xcms_object(data_container)
setGeneric(
  "is_xcms_object",
  function(object) standardGeneric("is_xcms_object")
)

#' @rdname is_xcms_object
setMethod(
  f = "is_xcms_object",
  signature = c("lcmsPlotDataContainer"),
  function(object) {
    is_xcms_data(object@data_obj)
    # inherits(object@data_obj, c("XCMSnExp", "MsExperiment"))
  }
)
