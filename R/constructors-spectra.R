#' Create a spectrum of the closest scan to the specified RT.
#'
#' @param raw_data An mzR object.
#' @param rt The RT to consider.
#' @param ms_level The MS level of the scans.
create_spectrum_from_closest_scan_to_rt <- function(raw_data, rt, ms_level) {
  hdr <- mzR::header(raw_data)

  ms_level_indices <- which(hdr$msLevel == ms_level)

  if (length(ms_level_indices) == 0) {
    stop("No scans found with the specified MS level.")
  }

  rt_diffs <- abs(hdr$retentionTime[ms_level_indices] - rt)

  closest_index <- ms_level_indices[which.min(rt_diffs)]
  closest_scan_id <- hdr[closest_index,]$seqNum

  spectrum_data <- mzR::peaks(raw_data, closest_scan_id)

  spectrum_df <- data.frame(
    mz = spectrum_data[, 1],
    intensity = spectrum_data[, 2],
    rt = hdr$retentionTime[closest_index]
  )

  return(spectrum_df)
}

#' Create a spectrum of the specified scan.
#'
#' @param raw_data An mzR object.
#' @param sample_metadata The sample metadata.
#' @param scan_index The scan index.
create_spectrum_from_scan_index <- function(
  raw_data,
  sample_metadata,
  scan_index
) {
  hdr <- mzR::header(raw_data)
  if (is.character(scan_index)) {
    sidx <- sample_metadata[[scan_index]]
  } else {
    sidx <- scan_index
  }
  spectrum_data <- mzR::peaks(raw_data, sidx)
  spectrum_df <- data.frame(
    mz = spectrum_data[, 1],
    intensity = spectrum_data[, 2],
    rt = hdr$retentionTime[sidx]
  )
  return(spectrum_df)
}

#' Create spectra for a single sample.
#'
#' @param raw_obj An mzR object.
#' @param detected_peaks The detected peaks to consider.
#' Only for modes 'closest_apex' and 'across_peak'.
#' @param sample_metadata The sample's metadata.
#' @param options The plot object's options.
#' @param rt_range The RT range to apply to the detected peaks.
#' @returns A data frame representing spectra
#' with columns 'mz', 'intensity', and 'rt'.
create_spectra_for_sample <- function(
  raw_obj,
  detected_peaks,
  sample_metadata,
  options,
  rt_range = NULL
) {
  spectra <- NULL
  sid <- sample_metadata$sample_id

  if (!is.null(options$spectra$scan_index)) {
    spectra <- create_spectrum_from_scan_index(
      raw_obj,
      sample_metadata,
      options$spectra$scan_index
    )
  } else if (options$spectra$mode == "closest" & !is.null(options$spectra$rt)) {
    # TODO: check if RT is outside range
    # If different RTs are supplied depending on the sample ID
    if (!is.null(names(options$spectra$rt))) {
      rt_to_consider <- unname(options$spectra$rt[sample_metadata$sample_id])
    } else {
      rt_to_consider <- options$spectra$rt
    }

    spectra <- create_spectrum_from_closest_scan_to_rt(
      raw_obj,
      rt = rt_to_consider,
      ms_level = options$spectra$ms_level
    )
  } else if (options$spectra$mode == "closest_apex") {
    spectra <- detected_peaks %>%
      filter(.data$sample_id == sid) %>%
      { if (!is.null(rt_range)) filter(., .data$rt >= rt_range[1], .data$rt <= rt_range[2]) else . } %>%
      pull(.data$rt) %>%
      lapply(function (rt) create_spectrum_from_closest_scan_to_rt(
        raw_obj,
        rt = rt,
        ms_level = options$spectra$ms_level)) %>%
      do.call(rbind, .)
  } else if (options$spectra$mode == "across_peak") {
    spectra <- detected_peaks %>%
      filter(.data$sample_id == sid) %>%
      { if (!is.null(rt_range)) filter(., .data$rt >= rt_range[1], .data$rt <= rt_range[2]) else . } %>%
      rowwise() %>%
      mutate(intervals = list(seq(.data$rtmin, .data$rtmax, by = options$spectra$interval))) %>%
      tidyr::unnest(.data$intervals) %>%
      mutate(rt_interval = .data$intervals) %>%
      select(.data$sample_id, .data$rt, .data$rtmin, .data$rtmax, .data$rt_interval) %>%
      pull(.data$rt_interval) %>%
      lapply(function (rt) create_spectrum_from_closest_scan_to_rt(
        raw_obj,
        rt = rt,
        ms_level = options$spectra$ms_level)) %>%
      do.call(rbind, .)
  }

  return(spectra)
}
