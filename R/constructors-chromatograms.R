create_bpc_tic <- function(raw_data, aggregation_fun) {
  hdr <- mzR::header(raw_data)
  ms1_header <- hdr[hdr$msLevel == 1, ]
  rt <- ms1_header$retentionTime

  if (aggregation_fun == "max") {
    bpi <- ms1_header$basePeakIntensity
  } else {
    bpi <- ms1_header$totIonCurrent
  }

  return(list(
    chromatograms = data.frame(rt = rt, intensity = bpi)
  ))
}

create_chromatogram <- function(raw_data, mz_range, rt_range) {
  hdr <- mzR::header(raw_data)
  scans_in_rt <- hdr[hdr$retentionTime >= rt_range[1] & hdr$retentionTime <= rt_range[2], ]

  chr <- data.frame()
  mass_traces <- data.frame()

  for (scan_id in scans_in_rt$seqNum) {
    spectrum <- mzR::peaks(raw_data, scan_id)
    rt <- scans_in_rt$retentionTime[scans_in_rt$seqNum == scan_id]
    in_mz_range <- spectrum[spectrum[, 1] >= mz_range[1] & spectrum[, 1] <= mz_range[2], , drop = FALSE]
    total_intensity <- sum(in_mz_range[, 2])

    if (nrow(in_mz_range) > 0) {
      chr <- rbind(chr, data.frame(rt = rt, intensity = total_intensity))
      mass_traces <- rbind(mass_traces, data.frame(rt = rt, mz = in_mz_range[, 1]))
    }
  }

  return(list(
    chromatograms = chr,
    mass_traces = mass_traces
  ))
}
