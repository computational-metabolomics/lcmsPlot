#' Get the metadata associated with the input object
#'
#' @param obj The data object
#' @param sample_id_column The column that should be associated with the sample ID
#' @param metadata The metadata data frame if not already included in the data object
#' @export
#' @examples
#' # Get metadata from sample paths
#' paths <- c("test01.mzML", "test02.mzML")
#' metadata <- get_metadata(paths, sample_id_column = NULL, metadata = NULL)
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
    MsExperiment::sampleData(obj) <- metadata
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
#' @examples
#' cdfs <- dir(
#'    system.file("cdf", package = "faahKO"),
#'    full.names = TRUE,
#'    recursive = TRUE)[c(1, 7)]
#' sample_names <- sub(basename(cdfs), pattern = ".CDF", replacement = "", fixed = TRUE)
#'
#' pd <- data.frame(sample_name = sample_names,
#'                  sample_group = c("KO", "WT"),
#'                  stringsAsFactors = FALSE)
#'
#' faahko <- readMsExperiment(spectraFiles = cdfs, sampleData = pd)
#'
#' cwp <- xcms::CentWaveParam(peakwidth = c(20, 80), noise = 5000, prefilter = c(6, 5000))
#' faahko <- xcms::findChromPeaks(faahko, param = cwp)
#'
#' detected_peaks <- get_detected_peaks(faahko)
get_detected_peaks <- function(obj) {
  UseMethod("get_detected_peaks")
}

.get_detected_peaks_xcms <- function(obj) {
  if (xcms::hasChromPeaks(obj)) {
    as.data.frame(xcms::chromPeaks(obj)) %>%
      dplyr::rename(sample_index = sample)
  } else {
    NULL
  }
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
#' @examples
#' cdfs <- dir(
#'    system.file("cdf", package = "faahKO"),
#'    full.names = TRUE,
#'    recursive = TRUE)[c(1, 7)]
#' sample_names <- sub(basename(cdfs), pattern = ".CDF", replacement = "", fixed = TRUE)
#'
#' pd <- data.frame(sample_name = sample_names,
#'                  sample_group = c("KO", "WT"),
#'                  stringsAsFactors = FALSE)
#'
#' faahko <- readMsExperiment(spectraFiles = cdfs, sampleData = pd)
#'
#' cwp <- xcms::CentWaveParam(peakwidth = c(20, 80), noise = 5000, prefilter = c(6, 5000))
#' faahko <- xcms::findChromPeaks(faahko, param = cwp)
#'
#' pdp <- xcms::PeakDensityParam(
#'     sampleGroups = pd$sample_group,
#'     minFraction = 0.4,
#'     bw = 30)
#' faahko <- xcms::groupChromPeaks(faahko, param = pdp)
#'
#' grouped_peaks <- get_grouped_peaks(faahko)
get_grouped_peaks <- function(obj) {
  UseMethod("get_grouped_peaks")
}

#' @rdname get_grouped_peaks
#' @export
get_grouped_peaks.default <- function(obj) {
  return(NULL)
}

.get_grouped_peaks_xcms <- function(obj) {
  if (xcms::hasFeatures(obj)) {
    as.data.frame(xcms::featureDefinitions(obj)) %>%
      rename(mz = mzmed, rt = rtmed) %>%
      mutate(name = xcms_utils$group_names(obj)) %>%
      xcms_utils$format_feature_identifiers(num_digits_rt = 0, num_digits_mz = 4)
  } else {
    NULL
  }
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
