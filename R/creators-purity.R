#' Populate the `purity_scores` slot from a `purityA` object
#'
#' Extracts per-scan precursor ion purity scores from `purityA@@puritydf` and
#' stores them in the `purity_scores` slot of the data container.
#'
#' @param obj An instance of class `lcmsPlotDataContainer` whose `data_obj`
#'   slot is a `purityA` object.
#' @param options A `list` of plot options.
#' @return A modified `lcmsPlotDataContainer` with the `purity_scores` slot
#'   populated.
#' @keywords internal
create_purity_scores <- function(obj, options) {
    pa <- obj@data_obj

    pdf <- as_tibble(pa@puritydf) |>
        dplyr::rename(
            rt = "precursorRT",
            in_purity = "inPurity",
            precursor_mz = "precursorMZ"
        ) |>
        dplyr::mutate(
            metadata_index = as.numeric(.data$fileid),
            feature_metadata_id = 0
        ) |>
        dplyr::select(
            "rt", "in_purity", "precursor_mz",
            "metadata_index", "feature_metadata_id"
        )

    sample_ids <- options$purity_scores_sample_ids
    if (!is.null(sample_ids)) {
        matching_idx <- obj@metadata |>
            dplyr::filter(.data$sample_id %in% sample_ids) |>
            dplyr::pull("sample_index")
        pdf <- pdf |> dplyr::filter(.data$metadata_index %in% matching_idx)
    }

    if (nrow(obj@chromatograms) > 0) {
        rt_min <- min(obj@chromatograms$rt, na.rm = TRUE)
        rt_max <- max(obj@chromatograms$rt, na.rm = TRUE)
        pdf <- pdf |> dplyr::filter(.data$rt >= rt_min & .data$rt <= rt_max)
    }

    obj@purity_scores <- pdf
    validObject(obj)
    return(obj)
}

#' Populate the `spectra` slot with an MS1 isolation window spectrum
#'
#' Extracts the MS1 survey scan nearest to a fragmentation event from
#' `purityA@@puritydf`, reads the spectrum from the raw file, and stores it
#' in the `spectra` slot. Isolation-window boundaries and the `inPurity` score
#' are stored in `feature_metadata` for downstream annotation by
#' `plot_isolation_window()`.
#'
#' @param obj An instance of class `lcmsPlotDataContainer` whose `data_obj`
#'   is a `purityA` object.
#' @param options A `list` of plot options (reads `options$isolation_window`).
#' @param pid An integer `pid` from `puritydf` selecting which fragmentation
#'   event(s) to visualise. `NULL` selects all events.
#' @param sample_id A character sample ID restricting which file to read.
#'   Ignored when `pid` is provided. `NULL` selects all samples.
#' @return A modified `lcmsPlotDataContainer` with `spectra` and
#'   `feature_metadata` slots populated.
#' @keywords internal
create_isolation_window <- function(obj, options, pid, sample_id) {
    pa <- obj@data_obj
    pdf <- as_tibble(pa@puritydf)

    if (!is.null(pid)) {
        selected <- pdf |> dplyr::filter(.data$pid %in% !!pid)
    } else if (!is.null(sample_id)) {
        meta <- obj@metadata
        idx <- meta |>
            dplyr::filter(.data$sample_id %in% !!sample_id) |>
            dplyr::pull("sample_index")
        selected <- pdf |> dplyr::filter(.data$fileid %in% idx)
    } else {
        selected <- pdf
    }

    if (nrow(selected) == 0) {
        stop("lp_isolation_window: no matching precursor records found in puritydf")
    }

    half_width <- options$isolation_window$half_width

    spectra_list <- vector("list", nrow(selected))
    feat_meta_list <- vector("list", nrow(selected))

    for (i in seq_len(nrow(selected))) {
        row <- selected[i, ]
        path <- pa@fileList[[row$fileid]]
        reader <- open_raw_reader(path)
        on.exit(ms_close(reader), add = TRUE)

        peaks_mat <- ms_peaks(reader, scans = row$precursorNearest)[[1L]]

        max_i <- max(peaks_mat[, 2], na.rm = TRUE)
        if (max_i == 0) {
            max_i <- 1
        }

        spectra_list[[i]] <- tibble(
            mz = peaks_mat[, 1],
            intensity = peaks_mat[, 2] / max_i * 100,
            rt = row$precursorRT,
            metadata_index = as.numeric(row$fileid),
            feature_metadata_id = as.numeric(i),
            reference = FALSE
        )

        feat_meta_list[[i]] <- tibble(
            feature_metadata_id = as.numeric(i),
            metadata_index = as.numeric(row$fileid),
            in_purity = row$inPurity,
            precursor_mz = row$precursorMZ,
            isolation_window_half_width = half_width
        )
    }

    obj@spectra <- do.call(rbind, spectra_list)
    obj@feature_metadata <- do.call(rbind, feat_meta_list)

    validObject(obj)
    return(obj)
}
