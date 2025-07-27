#' Get the metadata associated with the input object
#'
#' @param obj The data object
#' @param sample_id_column The column that should be associated with the sample ID
#' @param metadata The metadata data frame if not already included in the data object
#' @export
get_metadata <- function(obj, sample_id_column, metadata) {
  UseMethod("get_metadata")
}

#' @rdname get_metadata
#' @export
get_metadata.character <- function(obj, sample_id_column, metadata) {
  if (is.null(metadata)) {
    data.frame(sample_path = obj) %>%
      mutate(
        sample_index = row_number(),
        sample_id = tools::file_path_sans_ext(basename(obj))
      )
  } else {
    metadata %>%
      mutate(
        sample_index = row_number(),
        sample_id = .data[[sample_id_column]],
        sample_path = obj
      )
  }
}

#' @rdname get_metadata
#' @export
get_metadata.XCMSnExp <- function(obj, sample_id_column, metadata) {
  if (!is.null(metadata)) {
    xcms::phenoData(obj) <- new("NAnnotatedDataFrame", metadata)
  }

  xcms::phenoData(obj)@data %>%
    mutate(
      sample_index = row_number(),
      sample_id = .data[[sample_id_column]],
      sample_path = xcms::fileNames(obj)
    )
}

#' @rdname get_metadata
#' @export
get_metadata.MsExperiment <- function(obj, sample_id_column, metadata) {
  if (!is.null(metadata)) {
    MsExperiment::sampleData(obj) <- metadata # TODO: test this
  }

  MsExperiment::sampleData(obj) %>%
    as.data.frame() %>%
    mutate(
      sample_index = row_number(),
      sample_id = .data[[sample_id_column]],
      sample_path = xcms::fileNames(obj)
    )
}

#' Get the detected peaks from the data object (e.g. XCMSnExp)
#'
#' @param obj The data object
#' @export
get_detected_peaks <- function(obj) {
  UseMethod("get_detected_peaks")
}

.get_detected_peaks_xcms <- function(obj) {
  as.data.frame(xcms::chromPeaks(obj)) %>%
    dplyr::rename(sample_index = .data$sample)
}

#' @rdname get_detected_peaks
#' @export
get_detected_peaks.character <- function(obj) {
  return(NULL)
}

#' @rdname get_detected_peaks
#' @export
get_detected_peaks.XCMSnExp <- function(obj) {
  .get_detected_peaks_xcms(obj)
}

#' @rdname get_detected_peaks
#' @export
get_detected_peaks.MsExperiment <- function(obj) {
  .get_detected_peaks_xcms(obj)
}

#' Get the grouped peaks across samples (features) from the data object
#'
#' @param obj The data object
#' @export
get_grouped_peaks <- function(obj) {
  UseMethod("get_grouped_peaks")
}

#' @rdname get_grouped_peaks
#' @export
get_grouped_peaks.character <- function(obj) {
  return(NULL)
}

.get_grouped_peaks_xcms = function(obj) {
  as.data.frame(xcms::featureDefinitions(obj)) %>%
    rename(mz = .data$mzmed, rt = .data$rtmed) %>%
    mutate(name = xcms_utils$group_names(obj)) %>%
    xcms_utils$format_feature_identifiers(num_digits_rt = 0, num_digits_mz = 4)
}

#' @rdname get_grouped_peaks
#' @export
get_grouped_peaks.XCMSnExp <- function(obj) {
  .get_grouped_peaks_xcms(obj)
}

#' @rdname get_grouped_peaks
#' @export
get_grouped_peaks.XcmsExperiment <- function(obj) {
  .get_grouped_peaks_xcms(obj)
}

#' @rdname get_grouped_peaks
#' @export
get_grouped_peaks.MsExperiment <- function(obj) {
  .get_grouped_peaks_xcms(obj)
}
