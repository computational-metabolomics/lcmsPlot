#' Create an instance of class `lcmsPlotDataContainer` from peak counts per
#' retention-time bin
#'
#' This function creates an `lcmsPlotDataContainer` object holding the number of
#' chromatographic peaks each sample yielded in each retention-time bin,
#' mirroring `xcms::plotChromPeakImage()`.
#'
#' Bins with no peaks are kept with a count of zero rather than dropped, so a
#' sample that stopped producing peaks part-way through a run shows up as an
#' empty stretch instead of quietly disappearing from the plot.
#'
#' @param obj An instance of class `lcmsPlotDataContainer`.
#' @param options A `list` representing the plot object's options.
#' @return An instance of class `lcmsPlotDataContainer` with the created peak
#' count image. The `peak_count_image` slot is a `tibble` with columns `rt` (the
#' bin centre), `n_peaks`, `metadata_index`, and `feature_metadata_id`.
#' @keywords internal
setGeneric(
    "create_peak_count_image",
    function(obj, options) standardGeneric("create_peak_count_image")
)

#' @rdname create_peak_count_image
setMethod(
    f = "create_peak_count_image",
    signature = c("lcmsPlotDataContainer", "list"),
    definition = function(obj, options) {
        opts <- options$peak_count_image

        peaks <- get_detected_peaks(obj@data_obj)
        if (is.null(peaks) || nrow(peaks) == 0) {
            stop(
                "create_peak_count_image: the data object reports no ",
                "chromatographic peaks."
            )
        }
        peaks <- as_tibble(peaks)

        metadata <- obj@metadata
        if (!is.null(opts$sample_ids)) {
            metadata <- metadata |>
                filter(.data$sample_id %in% opts$sample_ids)
        }
        peaks <- peaks |>
            filter(.data$sample_index %in% metadata$sample_index)

        # The bin grid spans the whole acquisition, not just the range in which
        # peaks were found, so a short run reads as empty bins at the edge.
        rt_range <- opts$rt_range
        if (is.null(rt_range)) {
            rt_range <- tryCatch({
                rts <- xcms::rtime(obj@data_obj)
                if (length(rts) > 0) range(rts, na.rm = TRUE) else NULL
            }, error = function(e) NULL)
        }
        if (is.null(rt_range)) {
            rt_range <- range(peaks$rt, na.rm = TRUE)
        }
        rt_range <- range(rt_range)

        bin_size <- opts$bin_size
        breaks <- seq(floor(rt_range[1]), ceiling(rt_range[2]), by = bin_size)
        if (breaks[length(breaks)] < rt_range[2]) {
            breaks <- c(breaks, breaks[length(breaks)] + bin_size)
        }
        centres <- breaks[-length(breaks)] + bin_size / 2

        # cut(right = TRUE, include.lowest = TRUE) reproduces hist()'s binning,
        # which is what the xcms method uses.
        in_range <- peaks |>
            filter(
                .data$rt >= breaks[1],
                .data$rt <= breaks[length(breaks)]
            ) |>
            mutate(
                bin = cut(
                    .data$rt,
                    breaks = breaks,
                    include.lowest = TRUE,
                    labels = FALSE
                )
            )

        counts <- in_range |>
            count(.data$sample_index, .data$bin, name = "n_peaks")

        full_grid <- as_tibble(expand.grid(
            bin = seq_along(centres),
            sample_index = metadata$sample_index
        ))

        obj@peak_count_image <- full_grid |>
            left_join(counts, by = c("sample_index", "bin")) |>
            mutate(
                rt = centres[.data$bin],
                n_peaks = tidyr::replace_na(.data$n_peaks, 0),
                metadata_index = .data$sample_index,
                feature_metadata_id = NA_real_
            ) |>
            select(
                "rt", "n_peaks", "metadata_index", "feature_metadata_id"
            )

        validObject(obj)

        return(obj)
    }
)
