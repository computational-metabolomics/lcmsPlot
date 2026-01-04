#' Create a spectrum of the closest scan to the specified RT
#'
#' @param raw_data An instance of class `mzR`.
#' @param rt A `numeric` value indicating the RT to consider.
#' @param ms_level A `numeric` value indicating the MS level of the scans.
#' @return A `data.frame` representing a spectrum with columns
#' `mz`, `intensity`, and `rt`.
#' @keywords internal
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

#' Create a spectrum of the specified scan
#'
#' @param raw_data An instance of class `mzR`.
#' @param sample_metadata A `data.frame` indicating the sample metadata.
#' @param scan_index A `numeric` value indicating the scan index.
#' @return A `data.frame` representing a spectrum with columns
#' `mz`, `intensity`, and `rt`.
#' @keywords internal
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

#' Create spectra for a single sample
#'
#' @param raw_obj An instance of class `mzR`.
#' @param detected_peaks A `data.frame` containing the detected peaks to
#' consider. Only applicable for modes `"closest_apex"` and `"across_peak"`.
#' @param sample_metadata A `data.frame` containing the sample's metadata.
#' @param options A `list` representing the plot object's options.
#' @param rt_range A `numeric` value indicating the RT range to apply to the
#' detected peaks.
#' @return A `data.frame` representing spectra with columns `mz`, `intensity`,
#' and `rt`.
#' @keywords internal
create_spectra_for_sample <- function(
    raw_obj,
    detected_peaks,
    sample_metadata,
    options,
    rt_range = NULL
) {
    spectra <- NULL
    sid <- sample_metadata$sample_id
    rt_spectra <- options$spectra$rt

    if (!is.null(options$spectra$scan_index)) {
        spectra <- create_spectrum_from_scan_index(
            raw_obj,
            sample_metadata,
            options$spectra$scan_index
        )
    } else if (options$spectra$mode == "closest" & !is.null(rt_spectra)) {
        # TODO: check if RT is outside range
        # If different RTs are supplied depending on the sample ID
        if (!is.null(names(rt_spectra))) {
            rt_to_consider <- unname(rt_spectra[sample_metadata$sample_id])
        } else {
            rt_to_consider <- rt_spectra
        }

        spectra <- create_spectrum_from_closest_scan_to_rt(
            raw_obj,
            rt = rt_to_consider,
            ms_level = options$spectra$ms_level
        )
    } else if (options$spectra$mode == "closest_apex") {
        peaks <- detected_peaks |>
            filter(.data$sample_id == sid)

        if (!is.null(rt_range)) {
            peaks <- peaks |>
                filter(.data$rt >= rt_range[1], .data$rt <= rt_range[2])
        }

        spectra <- peaks |>
            pull(.data$rt) |>
            lapply(function(rt)
                create_spectrum_from_closest_scan_to_rt(
                    raw_obj,
                    rt = rt,
                    ms_level = options$spectra$ms_level
                )
            ) |>
            (\(x) do.call(rbind, x))()
    } else if (options$spectra$mode == "across_peak") {
        spectra <- detected_peaks |>
            filter(.data$sample_id == sid) |>
            (\(x) {
                if (!is.null(rt_range))
                    filter(x, .data$rt >= rt_range[1], .data$rt <= rt_range[2])
                else
                    x
            })() |>
            rowwise() |>
            mutate(intervals = list(seq(
                .data$rtmin,
                .data$rtmax,
                by = options$spectra$interval))
            ) |>
            tidyr::unnest(.data$intervals) |>
            mutate(rt_interval = .data$intervals) |>
            select(
                .data$sample_id,
                .data$rt,
                .data$rtmin,
                .data$rtmax,
                .data$rt_interval
            ) |>
            pull(.data$rt_interval) |>
            lapply(function(rt)
                create_spectrum_from_closest_scan_to_rt(
                    raw_obj,
                    rt = rt,
                    ms_level = options$spectra$ms_level
                )
            ) |>
            (\(x) do.call(rbind, x))()
    }

    return(spectra)
}
