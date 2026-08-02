#' Overlay chromatographic peak boundaries on an rt / m/z panel
#'
#' Draws one rectangle per detected chromatographic peak, spanning
#' `rtmin`-`rtmax` by `mzmin`-`mzmax`. This is the primitive behind
#' `xcms::plotChromPeaks()` and `xcms::plot(type = "XIC")`: the former draws the
#' rectangles on an empty frame, the latter over the ion map.
#'
#' The host renderer passes its own axis assignment, so the overlay follows when
#' the map is flipped to m/z on x.
#'
#' @param p A `ggplot` object whose axes are retention time and m/z.
#' @param detected_peaks A `data.frame` of detected peaks, requiring columns
#' `rtmin`, `rtmax`, `mzmin`, `mzmax`, and `sample_id`.
#' @param options A `list` of plot options (reads `options$chrom_peak_rects`).
#' @param x_dim A `character` naming the dimension on the x axis, `"rt"`
#' (default) or `"mz"`. The other axis carries the remaining dimension.
#' @param rt_range An optional length-2 `numeric` clamping the rectangles to the
#' host panel's retention-time window.
#' @param mz_range An optional length-2 `numeric` clamping the rectangles to the
#' host panel's m/z window.
#' @return The modified `ggplot` object.
#' @keywords internal
chrom_peak_rects_layer <- function(
    p,
    detected_peaks,
    options,
    x_dim = "rt",
    rt_range = NULL,
    mz_range = NULL
) {
    if (is.null(detected_peaks) || nrow(detected_peaks) == 0) {
        return(p)
    }

    required <- c("rtmin", "rtmax", "mzmin", "mzmax")
    missing_cols <- setdiff(required, names(detected_peaks))
    if (length(missing_cols) > 0) {
        stop(
            "lp_chrom_peak_rects: the detected peaks are missing the column(s) ",
            paste(missing_cols, collapse = ", "),
            ", so the peak boundaries cannot be drawn."
        )
    }

    opts <- options$chrom_peak_rects
    peaks <- detected_peaks

    if (!is.null(opts$sample_ids) && "sample_id" %in% names(peaks)) {
        peaks <- peaks[peaks$sample_id %in% opts$sample_ids, , drop = FALSE]
    }

    # Keep peaks that overlap the host window and clamp them to it, so the
    # overlay never widens the panel it is decorating.
    if (length(rt_range) == 2L) {
        rt_range <- range(rt_range)
        peaks <- peaks[
            peaks$rtmax >= rt_range[1] & peaks$rtmin <= rt_range[2], ,
            drop = FALSE]
        peaks$rtmin <- pmax(peaks$rtmin, rt_range[1])
        peaks$rtmax <- pmin(peaks$rtmax, rt_range[2])
    }
    if (length(mz_range) == 2L) {
        mz_range <- range(mz_range)
        peaks <- peaks[
            peaks$mzmax >= mz_range[1] & peaks$mzmin <= mz_range[2], ,
            drop = FALSE]
        peaks$mzmin <- pmax(peaks$mzmin, mz_range[1])
        peaks$mzmax <- pmin(peaks$mzmax, mz_range[2])
    }

    if (nrow(peaks) == 0) {
        return(p)
    }

    # On m/z-binned data (the faahKO CDFs, for instance) mzmin == mzmax and a
    # plain geom_rect draws nothing at all. Give every box a minimum height so a
    # degenerate peak still reads as a segment.
    span <- diff(range(c(peaks$mzmin, peaks$mzmax), na.rm = TRUE))
    min_height <- if (is.finite(span) && span > 0) span * 0.005 else 0.01
    pad <- pmax(0, (min_height - (peaks$mzmax - peaks$mzmin)) / 2)

    peaks$mz_lo <- peaks$mzmin - pad
    peaks$mz_hi <- peaks$mzmax + pad

    if (x_dim == "rt") {
        peaks$x_lo <- peaks$rtmin
        peaks$x_hi <- peaks$rtmax
        peaks$y_lo <- peaks$mz_lo
        peaks$y_hi <- peaks$mz_hi
    } else {
        peaks$x_lo <- peaks$mz_lo
        peaks$x_hi <- peaks$mz_hi
        peaks$y_lo <- peaks$rtmin
        peaks$y_hi <- peaks$rtmax
    }

    p + geom_rect(
        data = peaks,
        aes(
            xmin = .data$x_lo,
            xmax = .data$x_hi,
            ymin = .data$y_lo,
            ymax = .data$y_hi
        ),
        colour = opts$border,
        fill = opts$fill,
        alpha = opts$alpha,
        linewidth = opts$linewidth,
        inherit.aes = FALSE
    )
}
