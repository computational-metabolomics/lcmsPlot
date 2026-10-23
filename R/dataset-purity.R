#' Create the purity scores of a `purityA` source
#'
#' Extracts per-scan precursor ion purity scores from `purityA@@puritydf`,
#' restricted to the retention time range of the chromatograms when there are
#' any.
#'
#' @param obj An instance of class `lcmsPlotDataContainer` whose source is a
#'   `purityA` object.
#' @param options A `list` of plot options.
#' @return A modified `lcmsPlotDataContainer` with `summaries$purity_scores`
#'   set: a `tibble` with columns `rt`, `in_purity`, `precursor_mz`, and
#'   `sample_index`.
#' @keywords internal
create_purity_scores <- function(obj, options) {
    pa <- obj@source

    pdf <- as_tibble(pa@puritydf) |>
        dplyr::transmute(
            rt = .data$precursorRT,
            in_purity = .data$inPurity,
            precursor_mz = .data$precursorMZ,
            sample_index = as.integer(.data$fileid)
        )

    sample_ids <- options$purity_scores_sample_ids
    if (!is.null(sample_ids)) {
        samples <- .select_samples(obj, sample_ids)
        pdf <- pdf |>
            dplyr::filter(.data$sample_index %in% samples$sample_index)
    }

    rt_range <- .chromatogram_rt_range(obj@chromatograms)
    if (!is.null(rt_range)) {
        pdf <- pdf |>
            dplyr::filter(.data$rt >= rt_range[1] & .data$rt <= rt_range[2])
    }

    .set_summary(obj, "purity_scores", pdf)
}

# The retention time range covered by a set of chromatograms, or NULL when
# there are none.
.chromatogram_rt_range <- function(chromatograms) {
    if (length(chromatograms) == 0) {
        return(NULL)
    }
    rts <- unlist(
        lapply(Chromatograms::peaksData(chromatograms), `[[`, "rtime"),
        use.names = FALSE)
    if (length(rts) == 0) NULL else range(rts, na.rm = TRUE)
}

#' Populate the `spectra` slot with an MS1 isolation window spectrum
#'
#' Takes the MS1 survey scan nearest to each selected fragmentation event of
#' `purityA@@puritydf` and stores it in the `spectra` slot. The isolation
#' window (`precursor_mz`, `isolation_window_half_width`) and the `in_purity`
#' score are stored as spectra variables for `plot_isolation_window()`.
#'
#' @param obj An instance of class `lcmsPlotDataContainer` whose source is a
#'   `purityA` object.
#' @param options A `list` of plot options (reads `options$isolation_window`).
#' @param pid An integer `pid` from `puritydf` selecting which fragmentation
#'   event(s) to visualise. `NULL` selects all events.
#' @param sample_id A character sample ID restricting which file to read.
#'   Ignored when `pid` is provided. `NULL` selects all samples.
#' @return A modified `lcmsPlotDataContainer` with its `spectra` slot set.
#' @keywords internal
create_isolation_window <- function(obj, options, pid, sample_id) {
    pa <- obj@source
    pdf <- as_tibble(pa@puritydf)

    if (!is.null(pid)) {
        selected <- pdf |> dplyr::filter(.data$pid %in% !!pid)
    } else if (!is.null(sample_id)) {
        samples <- .select_samples(obj, sample_id)
        selected <- pdf |>
            dplyr::filter(.data$fileid %in% samples$sample_index)
    } else {
        selected <- pdf
    }

    if (nrow(selected) == 0) {
        stop("lp_isolation_window: no matching precursor records found in puritydf")
    }

    obj <- ensure_spectra(obj, options$parallel_param)
    targets <- tibble(
        sample_index = as.integer(selected$fileid),
        feature_id = as.character(seq_len(nrow(selected))),
        scan_index = as.integer(selected$precursorNearest)
    )

    sps <- select_target_spectra(obj, targets, ms_level = 1L)
    Spectra::rtime(sps) <- selected$precursorRT
    sps$in_purity <- selected$inPurity
    sps$precursor_mz <- selected$precursorMZ
    sps$isolation_window_half_width <- options$isolation_window$half_width

    obj@spectra <- scale_spectra(sps)
    obj
}
