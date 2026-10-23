#' Bring detected peaks to the common peak schema
#'
#' @param peaks A `data.frame` of detected peaks with at least `mz`, `rt`, and
#' `sample_index` columns.
#' @return A `tibble` with `mz`, `mzmin`, `mzmax`, `rt`, `rtmin`, `rtmax`,
#' `into`, `maxo`, and an integer `sample_index`, followed by any other column
#' of `peaks`. Missing boundary and intensity columns are filled with `NA`.
#' @keywords internal
normalise_peaks <- function(peaks) {
    if (is.null(peaks) || nrow(peaks) == 0) {
        return(tibble())
    }

    peaks <- as_tibble(peaks)
    for (col in c("mzmin", "mzmax", "rtmin", "rtmax", "into", "maxo")) {
        if (!col %in% colnames(peaks)) {
            peaks[[col]] <- NA_real_
        }
    }
    peaks$sample_index <- as.integer(peaks$sample_index)

    peaks |>
        relocate(all_of(c(
            "mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "into", "maxo",
            "sample_index")))
}

#' Keep the detected peaks that fall inside extraction windows
#'
#' A non-equi join of the peaks onto the windows of the same sample. Each
#' matched peak takes the `feature_id` of its window; a peak inside several
#' windows is kept once per window.
#'
#' @param peaks A `data.frame` of detected peaks, or `NULL`.
#' @param windows A `data.frame` of extraction windows (see
#' [extract_chromatograms()]).
#' @return A `data.frame` of the peaks inside the windows.
#' @keywords internal
peaks_in_windows <- function(peaks, windows) {
    if (is.null(peaks) || nrow(peaks) == 0 || nrow(windows) == 0) {
        return(tibble())
    }

    window_keys <- windows |>
        select("sample_index", "feature_id", "mzMin", "mzMax", "rtMin", "rtMax")

    peaks |>
        select(-any_of("feature_id")) |>
        mutate(sample_index = as.integer(.data$sample_index)) |>
        inner_join(
            window_keys |>
                mutate(sample_index = as.integer(.data$sample_index)),
            by = join_by(
                sample_index,
                between(mz, mzMin, mzMax),
                between(rt, rtMin, rtMax)),
            relationship = "many-to-many") |>
        select(-c("mzMin", "mzMax", "rtMin", "rtMax")) |>
        normalise_peaks()
}

# The detected peaks of the source, for the given samples.
.source_peaks <- function(container, samples = NULL) {
    peaks <- get_detected_peaks(container@source)
    if (is.null(peaks) || nrow(peaks) == 0) {
        return(tibble())
    }
    if (!is.null(samples)) {
        peaks <- peaks |> filter(.data$sample_index %in% samples$sample_index)
    }
    normalise_peaks(peaks)
}
