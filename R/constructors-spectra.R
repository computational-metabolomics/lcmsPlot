create_spectrum_from_closest_scan_to_rt <- function(raw_data, rt_range, rt, ms_level) {
  hdr <- mzR::header(raw_data)

  ms_level_indices <- which(hdr$msLevel == ms_level)

  if (length(ms_level_indices) == 0) {
    stop("No scans found with the specified MS level.")
  }

  rt_diffs <- abs(hdr$retentionTime[ms_level_indices] - rt)

  closest_index <- ms_level_indices[which.min(rt_diffs)]

  spectrum_data <- mzR::peaks(raw_data, closest_index)

  spectrum_df <- data.frame(
    mz = spectrum_data[, 1],
    intensity = spectrum_data[, 2],
    rt = hdr$retentionTime[closest_index]
  )

  return(spectrum_df)
}

create_spectrum_from_scan_index <- function(raw_data, sample_metadata, scan_index) {
  hdr <- mzR::header(raw_data)
  sidx <- ifelse(is.character(scan_index), sample_metadata[[scan_index]], scan_index)
  spectrum_data <- mzR::peaks(raw_data, sidx)
  spectrum_df <- data.frame(
    mz = spectrum_data[, 1],
    intensity = spectrum_data[, 2],
    rt = hdr$retentionTime[sidx]
  )
  return(spectrum_df)
}

create_spectra_for_sample <- function(raw_obj, detected_peaks, sample_metadata, options, rt_range = NULL) {
  spectra <- NULL
  sid <- sample_metadata$sample_id

  if (!is.null(options$spectra$scan_index)) {
    spectra <- create_spectrum_from_scan_index(raw_obj, sample_metadata, options$spectra$scan_index)
  } else if (options$spectra$mode == "closest" & !is.null(options$spectra$rt)) {
    # TODO: check if RT is outside range
    spectra <- create_spectrum_from_closest_scan_to_rt(
      raw_obj,
      rt = options$spectra$rt,
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
