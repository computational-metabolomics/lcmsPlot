#' Create a base peak or total ion current chromatogram
#'
#' @param raw_data An instance of class `MsRawReader`.
#' @param aggregation_fun A `character` value indicating the aggregation method.
#' One of `"max"` (base peak chromatogram), `"sum"` (total ion current), or
#' `"mean"` (averaged ion chromatogram).
#' @param rt_adjusted A `numeric` vector representing the adjusted RT values.
#' If `NULL` it will use the raw RT values.
#' @return A `list` with one `tibble` containing the chromatograms with
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

    intensity <- switch(
        aggregation_fun,
        max = ms1_header$basePeakIntensity,
        sum = ms1_header$totIonCurrent,
        mean = {
            mean_int <- ms1_header$meanIntensity
            if (is.null(mean_int) || all(is.na(mean_int))) {
                stop(
                    "aggregation_fun = \"mean\" is not available for ",
                    class(raw_data), ": the backend does not report a ",
                    "per-scan peak count. Use \"max\" or \"sum\" instead."
                )
            }
            mean_int
        },
        stop(
            "Unknown aggregation_fun: '", aggregation_fun, "'. ",
            "Expected 'max', 'sum', or 'mean'."
        )
    )

    return(list(
        chromatograms = tibble(rt = rt, intensity = intensity)
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
#' @param adjusted_rt A `tibble` containing the raw and adjusted RTs.
#' @param include_mass_traces A `logical` indicating whether to build the mass
#' traces. Callers that discard them should pass `FALSE`, which on the `rawrr`
#' backend avoids reading the scan header and spectra.
#' @return A `list` with two data frames (`chromatograms` and `mass_traces`)
#' containing the chromatograms with columns `rt` and `intensity`
#' and mass traces with columns `rt` and `mz`. `mass_traces` is empty when
#' `include_mass_traces` is `FALSE`.
#' @keywords internal
create_chromatogram <- function(
    raw_data,
    mz_range,
    rt_range,
    ms_level = 1,
    fill_gaps = FALSE,
    adjusted_rt = NULL,
    include_mass_traces = TRUE
) {
    # Guard against an undefined extraction window (e.g. a compound with no RT or
    # m/z): an NA range would otherwise corrupt the scan filtering below.
    if (any(is.na(mz_range)) || any(is.na(rt_range))) {
        return(list(
            chromatograms = tibble(rt = numeric(), intensity = numeric()),
            mass_traces = tibble(rt = numeric(), mz = numeric())
        ))
    }

    if (is(raw_data, "RawrrReader")) {
        xic <- .mz_range_to_xic(mz_range)
        chr <- ms_chromatogram(raw_data, xic$mz, xic$ppm, rt_range)[[1]]

        # Only enable mass traces when there are less than 100 scans.
        if (include_mass_traces && nrow(chr) <= 100) {
            hdr <- ms_header(raw_data)
            scans_in_rt <- hdr[
                hdr$retentionTime >= rt_range[1] &
                    hdr$retentionTime <= rt_range[2],
            ]
            spectra <- ms_peaks(raw_data, scans_in_rt$seqNum)
            mass_traces <- tibble()

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
                        tibble(rt = rt, mz = in_mz_range[, 1])
                    )
                }
            }
        } else {
            mass_traces <- tibble(
                rt = numeric(),
                mz = numeric()
            )
        }
    } else if (is(raw_data, "MzrReader") || is(raw_data, "XcmsRawReader")) {
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

        chr <- tibble(rt = numeric(), intensity = numeric())
        mass_traces <- tibble()

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
                chr <- rbind(chr, tibble(rt = rt, intensity = total_intensity))
                if (include_mass_traces) {
                    mass_trace <- tibble(rt = rt, mz = in_mz_range[, 1])
                    mass_traces <- rbind(mass_traces, mass_trace)
                }
            } else if (fill_gaps) {
                chr <- rbind(chr, tibble(rt = rt, intensity = 0))
            }
        }

        if (!include_mass_traces) {
            mass_traces <- tibble(rt = numeric(), mz = numeric())
        }
    } else {
        stop("Input raw data is not of a supported type.")
    }

    return(list(
        chromatograms = chr,
        mass_traces = mass_traces
    ))
}

#' Create extracted ion chromatograms for several windows of one file
#'
#' Batched counterpart of [create_chromatogram()] for callers that only need
#' the chromatograms. On the `rawrr` backend all windows are extracted in a
#' single read of the `.raw` file; other backends extract each window in turn.
#' Mass traces are not built.
#'
#' @param raw_data An instance of class `MsRawReader`.
#' @param mz_ranges A two-column `matrix` with one m/z range per row.
#' @param rt_ranges A two-column `matrix` with one RT range per row.
#' @param fill_gaps A `logical` indicating whether to fill gaps
#' between scans with zeros.
#' @return A `list` with one `tibble` (columns `rt` and `intensity`) per row
#' of `mz_ranges`. Windows with an undefined range give an empty `tibble`.
#' @keywords internal
create_chromatograms_batch <- function(
    raw_data,
    mz_ranges,
    rt_ranges,
    fill_gaps = FALSE
) {
    n <- nrow(mz_ranges)
    result <- rep(
        list(tibble(rt = numeric(), intensity = numeric())), n)
    valid <- which(
        rowSums(is.na(mz_ranges)) == 0 & rowSums(is.na(rt_ranges)) == 0)

    if (length(valid) == 0) {
        return(result)
    }

    if (is(raw_data, "RawrrReader")) {
        xics <- lapply(valid, function(i) .mz_range_to_xic(mz_ranges[i, ]))
        result[valid] <- ms_chromatogram(
            raw_data,
            mz = vapply(xics, `[[`, numeric(1), "mz"),
            ppm = vapply(xics, `[[`, numeric(1), "ppm"),
            rt_range = rt_ranges[valid, , drop = FALSE]
        )
    } else {
        result[valid] <- lapply(valid, function(i) {
            create_chromatogram(
                raw_data,
                mz_range = mz_ranges[i, ],
                rt_range = rt_ranges[i, ],
                fill_gaps = fill_gaps,
                include_mass_traces = FALSE
            )$chromatograms
        })
    }

    result
}

# Convert an m/z range into the centre m/z and ppm tolerance that
# rawrr::readChromatogram() expects.
.mz_range_to_xic <- function(mz_range) {
    mz <- (mz_range[1] + mz_range[2]) / 2
    list(mz = mz, ppm = ((mz_range[2] - mz_range[1]) / mz) * 1e6)
}
