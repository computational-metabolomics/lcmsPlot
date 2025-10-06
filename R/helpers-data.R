#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL

#' Check whether an object is from XCMS.
#'
#' @param obj The input object.
#' @returns Whether the object is an XCMS one.
is_xcms_data <- function(obj) {
  inherits(obj, c("XCMSnExp", "MsExperiment"))
}

#' Check whether an object is from XCMS and processed (e.g., peak detection).
#'
#' @param obj The input object.
#' @returns Whether the object is an XCMS one.
is_xcms_processed_data <- function(obj) {
  inherits(obj, c("XCMSnExp", "XcmsExperiment"))
}

#' Get an XCMSnExp example object
#'
#' @param indices The faahKO package data samples indices
#' @param should_group_peaks Whether to group the detected peaks
#' @returns The XCMSnExp object
#' @export
#' @examples
#' get_XCMSnExp_object_example(indices = 1:5)
get_XCMSnExp_object_example <- function(
  indices = 1:3,
  should_group_peaks = FALSE
) {
  cdfs <- dir(system.file("cdf", package = "faahKO"), full.names = TRUE,
              recursive = TRUE)[indices]
  sample_names <- sub(
    basename(cdfs),
    pattern = ".CDF",
    replacement = "",
    fixed = TRUE
  )

  pd <- data.frame(
    sample_name = sample_names,
    sample_group = toupper(sub("[0-9]+", "", sample_names)),
    stringsAsFactors = FALSE
  )

  raw_data <- MSnbase::readMSData(
    files = cdfs,
    pdata = new("AnnotatedDataFrame", pd),
    mode = "onDisk",
    msLevel = 1
  )

  cwp <- xcms::CentWaveParam(
    peakwidth = c(20, 80),
    noise = 10000,
    prefilter = c(6, 10000)
  )
  xdata <- xcms::findChromPeaks(raw_data, param = cwp)

  if (should_group_peaks) {
    xdata <- xcms::adjustRtime(xdata, param = xcms::ObiwarpParam(binSize = 0.6))
    pdp <- xcms::PeakDensityParam(
      sampleGroups = pd$sample_group,
      minFraction = 1,
      bw = 30)
    xdata <- xcms::groupChromPeaks(xdata, param = pdp)
  }

  xdata
}

#' Merge two data frames by an index column.
#'
#' @param a The first data frame.
#' @param b The second data frame.
#' @param index_col The indexing column in data frame 'b'.
#' @returns The merged data frame.
merge_by_index <- function(a, b, index_col) {
  b_mod <- b %>%
    mutate(row_id = row_number())

  a %>%
    left_join(b_mod, by = setNames("row_id", index_col))
}

#' Remove NULL values from a list.
#'
#' @param lst The input list.
#' @returns The list without NULL values.
remove_null_elements <- function(lst) {
  return(lst[!sapply(lst, is.null)])
}

#' Get the m/z range for an input m/z and a ppm error.
#'
#' @param mz The m/z value to calculate the range for.
#' @param ppm The ppm tolerance for the range.
#' @returns A vector with two values indicating the m/z range.
get_mz_range <- function(mz, ppm = 5) {
  mzdev <- mz * (ppm / 1000000)
  return(c(mz - mzdev, mz + mzdev))
}

#' Get a feature's m/z and RT ranges given different feature specifications.
#'
#' @param feature The input feature which can be a vector or a data frame row.
#' The accepted columns combinations are:
#' - mz, rt (optional)
#' - mzmin, mzmax, rtmin (optional), rtmax (optional)
#' @param options The plot object's options.
#' @param full_rt_range The full RT range if an RT range is not given.
#' @returns A named list defining a feature with names: feature_id, mzr, rtr.
get_feature_data <- function(feature, options, full_rt_range) {
  # Helper: safely extract a value by name from vector or data.frame
  get_val <- function(x, name) {
    if (is.data.frame(x)) {
      if (name %in% names(x)) return(x[[name]][1])
    } else if (is.vector(x)) {
      if (name %in% names(x)) return(x[[name]])
    }
    return(NA)
  }

  mz     <- get_val(feature, "mz")
  rt     <- get_val(feature, "rt")
  mzmin  <- get_val(feature, "mzmin")
  mzmax  <- get_val(feature, "mzmax")
  rtmin  <- get_val(feature, "rtmin")
  rtmax  <- get_val(feature, "rtmax")

  # Case 1: Single mz (± tolerance) with optional rt
  if (!is.na(mz)) {
    mzr <- get_mz_range(mz, options$chromatograms$ppm)

    if (!is.na(rt)) {
      rtr <- c(rt - options$chromatograms$rt_tol,
               rt + options$chromatograms$rt_tol)
      feature_id <- sprintf("M%dT%d", round(mz), round(rt))
    } else {
      rtr <- full_rt_range
      feature_id <- sprintf("M%d", round(mz))
    }

    # Case 2: Explicit mz range, with optional rt range
  } else if (!is.na(mzmin) && !is.na(mzmax)) {
    mzr <- c(mzmin, mzmax)

    if (!is.na(rtmin) && !is.na(rtmax)) {
      rtr <- c(rtmin, rtmax)
      feature_id <- sprintf("M%dT%d",
                            round(mean(c(mzmin, mzmax))),
                            round(mean(c(rtmin, rtmax))))
    } else {
      rtr <- full_rt_range
      feature_id <- sprintf("M%d", round(mean(c(mzmin, mzmax))))
    }

  } else {
    stop("Unsupported feature format: must provide either 'mz' or ('mzmin' and 'mzmax').")
  }

  list(
    feature_id = feature_id,
    mzr = mzr,
    rtr = rtr
  )
}

#' Get the feature data (m/z and RT ranges) for a set of features.
#'
#' @param options The plot object's options. This contains the input features.
#' @param sample_metadata The sample's metadata.
#' @param grouped_peaks The grouped peaks to use as input features.
#' @param full_rt_range The full RT range if an RT range is not given.
#' @returns A list of feature data (see get_feature_data).
get_features <- function(
  options,
  sample_metadata,
  grouped_peaks = NULL,
  full_rt_range = NULL
) {
  if (is.null(grouped_peaks)) {
    input_features <- options$chromatograms$features
  } else {
    input_features <- grouped_peaks
  }

  if (is.data.frame(input_features) &&
      "sample_id" %in% colnames(input_features)) {
    feature_indices <- which(input_features$sample_id == sample_metadata$sample_id)
  } else {
    feature_indices <- seq_len(nrow(input_features))
  }

  features <- list()

  for (i in seq_len(length(feature_indices))) {
    feature_index <- feature_indices[i]
    feature <- input_features[feature_index, ]
    features[[i]] <- get_feature_data(feature, options, full_rt_range)
  }

  return(features)
}
