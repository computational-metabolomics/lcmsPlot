create_spectrum_from_closest_scan_to_rt <- function(raw_data, rt, ms_level) {
  hdr <- mzR::header(raw_data)

  ms_level_indices <- which(hdr$msLevel == ms_level)

  if (length(ms_level_indices) == 0) {
    # TODO: rethink this
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
