#' @export
get_metadata <- function(obj, sample_id_column, metadata) {
  UseMethod("get_metadata")
}

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

#' @export
get_metadata.XCMSnExp <- function(obj, sample_id_column, metadata) {
  if (!is.null(metadata)) {
    phenoData(obj) <- new("NAnnotatedDataFrame", metadata)
  }

  phenoData(obj)@data %>%
    mutate(
      sample_index = row_number(),
      sample_id = .data[[sample_id_column]],
      sample_path = fileNames(obj)
    )
}

#' @export
get_detected_peaks <- function(obj) {
  UseMethod("get_detected_peaks")
}

#' @export
get_detected_peaks.XCMSnExp <- function(obj) {
  as.data.frame(chromPeaks(obj))
}

#' @export
get_grouped_peaks <- function(obj) {
  UseMethod("get_grouped_peaks")
}

#' @export
get_grouped_peaks.XCMSnExp <- function(obj) {
  as.data.frame(featureDefinitions(obj)) %>%
    rename(mz = mzmed, rt = rtmed) %>%
    mutate(name = groupnames(obj)) %>%
    xcms_utils$format_feature_identifiers(num_digits_rt = 0, num_digits_mz = 4)
}
