#' Create a base peak or total ion current chromatogram
#'
#' @param raw_data An instance of class `MsRawReader`.
#' @param aggregation_fun A `function` indicating the aggregation method.
#' One of `"sum"` or `"max"`.
#' @param rt_adjusted A `numeric` vector representing the adjusted RT values.
#' If `NULL` it will use the raw RT values.
#' @return A `list` with one `data.frame` containing the chromatograms with
#' columns `rt` and `intensity`.
#' @keywords internal
create_bpc_tic <- function(raw_data, aggregation_fun, rt_adjusted = NULL) {
    hdr <- ms_header(raw_data)
    ms1_header <- hdr[hdr$msLevel == 1, ]

    if (is.null(rt_adjusted)) {
        rt <- ms1_header$retentionTime
    } else {
        rt <- rt_adjusted
    }

    if (aggregation_fun == "max") {
        bpi <- ms1_header$basePeakIntensity
    } else {
        bpi <- ms1_header$totIonCurrent
    }

    return(list(
        chromatograms = data.frame(rt = rt, intensity = bpi)
    ))
}

#' Create an extracted ion chromatogram
#'
#' @param raw_data An instance of class `MsRawReader`.
#' @param mz_range A `numeric` vector indicating the m/z range.
#' @param rt_range A `numeric` vector indicating the RT range.
#' @param ms_level A `numeric` value indicating the MS level
#' of the scans to consider.
#' @param fill_gaps A `logical` indicating whether to fill gaps
#' between scans with zeros.
#' @param adjusted_rt A `data.frame`
#' @return A `list` with two data frames (`chromatograms` and `mass_traces`)
#' containing the chromatograms with columns `rt` and `intensity`
#' and mass traces with columns `rt` and `mz`.
#' @keywords internal
create_chromatogram <- function(
    raw_data,
    mz_range,
    rt_range,
    ms_level = 1,
    fill_gaps = FALSE,
    adjusted_rt = NULL
) {
    if (is(raw_data, "RawrrReader")) {
        mz <- (mz_range[1] + mz_range[2]) / 2
        ppm <- ((mz_range[2] - mz_range[1]) / mz) * 1e6
        chr <- ms_chromatogram(raw_data, mz, ppm, rt_range)

        # Only enable mass traces when there are less than 100 scans.
        if (nrow(chr) <= 100) {
            hdr <- ms_header(raw_data)
            scans_in_rt <- hdr[
                hdr$retentionTime >= rt_range[1] &
                    hdr$retentionTime <= rt_range[2],
            ]
            spectra <- ms_peaks(raw_data, scans_in_rt$seqNum)
            mass_traces <- data.frame()

            for (i in seq_len(nrow(scans_in_rt))) {
                spectrum <- spectra[[i]]
                rt <- scans_in_rt[i, ]$retentionTime
                in_mz_range <- spectrum[
                    spectrum[, 1] >= mz_range[1] &
                        spectrum[, 1] <= mz_range[2], ,
                    drop = FALSE
                ]
                if (nrow(in_mz_range) > 0) {
                    mass_traces <- rbind(
                        mass_traces,
                        data.frame(rt = rt, mz = in_mz_range[, 1])
                    )
                }
            }
        } else {
            mass_traces <- data.frame(
                rt = numeric(),
                mz = numeric()
            )
        }
    } else if (is(raw_data, "MzrReader")) {
        hdr <- ms_header(raw_data)

        if (!is.null(adjusted_rt) && nrow(adjusted_rt) == nrow(hdr)) {
            hdr$retentionTime <- adjusted_rt |> pull(.data$adj_rt)
        }

        hdr <- hdr[hdr$msLevel == ms_level, ]
        scans_in_rt <- hdr[
            hdr$retentionTime >= rt_range[1] &
                hdr$retentionTime <= rt_range[2],
        ]
        spectra <- ms_peaks(raw_data, scans_in_rt$seqNum)

        chr <- data.frame()
        mass_traces <- data.frame()

        for (i in seq_len(nrow(scans_in_rt))) {
            spectrum <- spectra[[i]]
            rt <- scans_in_rt[i, ]$retentionTime
            in_mz_range <- spectrum[
                spectrum[, 1] >= mz_range[1] &
                    spectrum[, 1] <= mz_range[2], ,
                drop = FALSE
            ]
            total_intensity <- sum(in_mz_range[, 2])

            if (nrow(in_mz_range) > 0) {
                chr <- rbind(chr, data.frame(rt = rt, intensity = total_intensity))
                mass_trace <- data.frame(rt = rt, mz = in_mz_range[, 1])
                mass_traces <- rbind(mass_traces, mass_trace)
            } else if (fill_gaps) {
                chr <- rbind(chr, data.frame(rt = rt, intensity = 0))
            }
        }
    } else {
        stop("Input raw data is not of a supported type.")
    }


    return(list(
        chromatograms = chr,
        mass_traces = mass_traces
    ))
}
