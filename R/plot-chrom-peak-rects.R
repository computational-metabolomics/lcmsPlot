#' Draw chromatographic peak boundaries on an rt / m/z panel
#'
#' Draws one rectangle per detected chromatographic peak, spanning
#' `rtmin`-`rtmax` by `mzmin`-`mzmax`.
#'
#' The caller passes its own axis assignment, so the rectangles follow when the
#' host map is flipped to m/z on x. Used both by [plot_chrom_peak_rects()] for
#' the standalone panel and by the mass trace and intensity map renderers when
#' the layer decorates them.
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

#' Plot chromatographic peak boundaries on their own panel
#'
#' It is used when `lp_chrom_peak_rects()` is called without a host layer; with
#' `lp_mass_trace()` or `lp_intensity_map()` present the rectangles are drawn as
#' an overlay by those renderers instead.
#'
#' @param datasets A named `list` of data frames containing the primary datasets
#' to plot. The used key is `chrom_peak_rects`.
#' @param supporting_datasets A `list` of supporting data frames. Unused;
#' present for API consistency.
#' @param options A list of plot options, controlling units, faceting,
#' highlighting, and other visual parameters.
#' @param single A `logical` value that indicates whether it should treat
#' the plot as a single dataset variant, which can affect faceting
#' and layout behavior. Default is `FALSE`.
#' @return A `ggplot` object representing the peak boundary panel.
#' @keywords internal
plot_chrom_peak_rects <- function(
    datasets,
    supporting_datasets,
    options,
    single = FALSE
) {
    dataset <- datasets$chrom_peak_rects
    opts <- options$chrom_peak_rects

    if (!is.null(opts$sample_ids) && "sample_id" %in% names(dataset)) {
        dataset <- dataset[dataset$sample_id %in% opts$sample_ids, , drop = FALSE]
    }

    extra_layers <- list(
        legend_title(options),
        faceting(options, single),
        grid_layout(options, single)
    )
    extra_layers <- remove_null_elements(extra_layers)

    p <- ggplot(dataset, aes(x = .data$rt, y = .data$mz)) +
        labs(x = "RT (sec)", y = "m/z") +
        theme_minimal() +
        extra_layers

    if (is.null(options$facets$facets) &&
        "sample_id" %in% names(dataset) &&
        length(unique(dataset$sample_id)) > 1) {
        p <- p + facet_wrap(~ sample_id)
    }

    chrom_peak_rects_layer(
        p,
        dataset,
        options,
        rt_range = opts$rt_range,
        mz_range = opts$mz_range)
}
