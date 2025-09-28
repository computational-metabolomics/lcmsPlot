#' Create a base peak or total ion current chromatogram.
#' @param raw_data An mzR object.
#' @param aggregation_fun The aggregation function to create the chromatogram.
#' One of 'sum' or 'max'.
#' @param rt_adjusted The RT adjusted values to use instead of the raw ones.
#' @returns A list with one data frame containing the chromatograms.
create_bpc_tic <- function(raw_data, aggregation_fun, rt_adjusted = NULL) {
  hdr <- mzR::header(raw_data)
  ms1_header <- hdr[hdr$msLevel == 1, ]

  if (is.null(rt_adjusted)) {
    rt <- ms1_header$retentionTime
  } else {
    rt <- rt_adjusted
  }

  if (is(raw_data, "mzRnetCDF")) {
    scan_nums <- ms1_header$seqNum
    bpi <- numeric(length(scan_nums))

    batch_size <- 100
    num_scans <- length(scan_nums)
    batch_starts <- seq(1, num_scans, by = batch_size)

    for (start in batch_starts) {
      end <- min(start + batch_size - 1, num_scans)
      batch_scans <- scan_nums[start:end]

      peaks_list <- mzR::peaks(raw_data, scans = batch_scans)

      for (i in seq_along(batch_scans)) {
        pk <- peaks_list[[i]]
        if (aggregation_fun == "max") {
          bpi[start + i - 1] <- max(pk[, 2])
        } else {
          bpi[start + i - 1] <- sum(pk[, 2])
        }
      }
    }
  } else {
    if (aggregation_fun == "max") {
      bpi <- ms1_header$basePeakIntensity
    } else {
      bpi <- ms1_header$totIonCurrent
    }
  }

  return(list(
    chromatograms = data.frame(rt = rt, intensity = bpi)
  ))
}

#' Create an extracted ion chromatogram.
#'
#' @param raw_data An mzR object.
#' @param mz_range The m/z range of the chromatogram.
#' @param rt_range The RT range of the chromatogram.
#' @param ms_level The MS level of the scans to consider.
#' @param fill_gaps Whether to fill gaps between scans with zeros.
create_chromatogram <- function(
  raw_data,
  mz_range,
  rt_range,
  ms_level = 1,
  fill_gaps = FALSE
) {
  hdr <- mzR::header(raw_data)
  hdr <- hdr[hdr$msLevel == ms_level,]
  scans_in_rt <- hdr[
    hdr$retentionTime >= rt_range[1] &
      hdr$retentionTime <= rt_range[2], ]
  spectra <- mzR::peaks(raw_data, scans_in_rt$seqNum)

  chr <- data.frame()
  mass_traces <- data.frame()

  for (i in seq_len(nrow(scans_in_rt))) {
    spectrum <- spectra[[i]]
    rt <- scans_in_rt[i,]$retentionTime
    in_mz_range <- spectrum[
      spectrum[, 1] >= mz_range[1] &
        spectrum[, 1] <= mz_range[2], , drop = FALSE]
    total_intensity <- sum(in_mz_range[, 2])

    if (nrow(in_mz_range) > 0) {
      chr <- rbind(chr, data.frame(rt = rt, intensity = total_intensity))
      mass_trace <- data.frame(rt = rt, mz = in_mz_range[, 1])
      mass_traces <- rbind(mass_traces, mass_trace)
    } else if (fill_gaps) {
      chr <- rbind(chr, data.frame(rt = rt, intensity = 0))
    }
  }

  return(list(
    chromatograms = chr,
    mass_traces = mass_traces
  ))
}
