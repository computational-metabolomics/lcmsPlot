#' Plot chromatographic peak counts per retention-time bin
#'
#' `plot_peak_count_image()` draws retention-time bins on x and samples on y,
#' filled by how many chromatographic peaks each sample yielded in each bin. It
#' is the ggplot2 counterpart of `xcms::plotChromPeakImage()`, and is read as a
#' quality-control view: an empty stretch marks a sample that stopped producing
#' peaks, a pale column a retention-time region where detection collapsed.
#'
#' @param datasets A named `list` of data frames containing the primary datasets
#' to plot. The used key is `peak_count_image`.
#' @param supporting_datasets A `list` of supporting data frames. Unused;
#' present for API consistency.
#' @param options A list of plot options, controlling units, faceting,
#' highlighting, and other visual parameters.
#' @param single A `logical` value that indicates whether it should treat
#' the plot as a single dataset variant, which can affect faceting
#' and layout behavior. Default is `FALSE`.
#' @return A `ggplot` object representing the peak count image.
#' @keywords internal
plot_peak_count_image <- function(
    datasets,
    supporting_datasets,
    options,
    single = FALSE
) {
    dataset <- datasets$peak_count_image
    opts <- options$peak_count_image

    extra_layers <- list(
        legend_title(options),
        faceting(options, single),
        grid_layout(options, single)
    )
    extra_layers <- remove_null_elements(extra_layers)

    # Injection order, not alphabetical: a run drifting over time is only
    # readable when the samples stay in the order they were acquired.
    if (all(c("sample_id", "sample_index") %in% names(dataset))) {
        levels_in_order <- dataset |>
            distinct(.data$sample_id, .data$sample_index) |>
            arrange(.data$sample_index) |>
            pull("sample_id")
        dataset$sample_label <- factor(
            dataset$sample_id,
            levels = levels_in_order)
    } else {
        dataset$sample_label <- factor(dataset$metadata_index)
    }

    if (isTRUE(opts$log)) {
        # log2(0) is -Inf, which would swallow the whole colour scale; leave
        # empty bins blank instead. xcms passes the -Inf straight to image().
        fill_values <- log2(dataset$n_peaks)
        fill_values[!is.finite(fill_values)] <- NA_real_
        fill_label <- "log2(peaks)"
    } else {
        fill_values <- dataset$n_peaks
        fill_label <- "Peaks"
    }
    dataset$fill_value <- fill_values

    p <- ggplot(
        dataset,
        aes(
            x = .data$rt,
            y = .data$sample_label,
            fill = .data$fill_value)
    ) +
        geom_tile() +
        (if (!is.null(opts$fill_scale)) {
            opts$fill_scale
        } else {
            scale_fill_viridis_c()
        }) +
        labs(x = "RT (sec)", y = NULL, fill = fill_label) +
        theme_minimal() +
        extra_layers

    return(p)
}
