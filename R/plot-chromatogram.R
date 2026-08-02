#' Plot chromatograms from one or more datasets
#'
#' `plot_chromatogram()` generates a ggplot2 chromatogram plot from processed
#' datasets. It can handle multiple datasets, optionally highlight detected
#' peaks or apices, and apply faceting or grid layouts based on plot options.
#'
#' @param datasets A named `list` of data frames containing the primary datasets
#' to plot. Typically includes `chromatograms` and optionally `spectra`.
#' @param supporting_datasets A `list` of supporting data frames, such as
#' detected peaks, used for highlighting features.
#' @param options A list of plot options, controlling units, faceting,
#' highlighting, and other visual parameters.
#' @param single A `logical` value that indicates whether it should treat
#' the plot as a single dataset variant, which can affect faceting
#' and layout behavior. Default is `FALSE`.
#' @return A `ggplot` object representing the chromatogram plot.
#' @keywords internal
plot_chromatogram <- function(
    datasets,
    supporting_datasets,
    options,
    single = FALSE
) {
    dataset <- datasets$chromatograms

    grouping_vars <- get_grouping_variables(options)

    # One line per sample, split further by the arrangement column. Colouring by
    # anything other than sample_id (feature_id being the obvious case) would
    # otherwise draw every series in a panel as a single zig-zagging line.
    series_cols <- intersect(
        unique(c("sample_id", options$arrangement$group_by)),
        names(dataset))
    dataset$series_id <- if (length(series_cols) > 0) {
        do.call(
            paste,
            c(unname(as.list(dataset[series_cols])), list(sep = "\r")))
    } else {
        "1"
    }

    dataset <- dataset |>
        group_by(across(all_of(grouping_vars))) |>
        mutate(
            rt_plot = case_when(
                options$chromatograms$rt_unit == "minute" ~ rt / 60,
                TRUE ~ rt),
            intensity_plot = case_when(
                options$chromatograms$intensity_unit == "relative" ~
                    (intensity / max(intensity)) * 100,
                TRUE ~ intensity)
        ) |>
        ungroup()

    transform_fun <- options$chromatograms$transform
    if (!is.function(transform_fun)) transform_fun <- identity

    stacked <- options$chromatograms$stacked
    if (!is.numeric(stacked) || length(stacked) == 0L) stacked <- 0
    stacked <- stacked[1L]

    # xcms derives the stacking band from transform(range(c(intensity, 0))) and
    # collapses the -Inf that log-like transforms produce at zero.
    y_limits <- transform_fun(c(0, max(dataset$intensity_plot, na.rm = TRUE)))
    y_limits[!is.finite(y_limits)] <- 0

    dataset$intensity_plot <- transform_fun(dataset$intensity_plot)

    # The series overlaid inside one panel are the ones the arrangement colours;
    # stacking them per (sample, feature) pair instead would make panels
    # incomparable.
    stack_col <- options$arrangement$group_by
    if (is.null(stack_col) || !stack_col %in% names(dataset)) {
        stack_col <- if ("sample_id" %in% names(dataset)) "sample_id" else NULL
    }

    offsets <- NULL
    if (stacked != 0 && !is.null(stack_col)) {
        offsets <- stacked_offsets(dataset, stack_col, stacked, y_limits)
        dataset$intensity_plot <- dataset$intensity_plot +
            series_offset(dataset, offsets, stack_col)
    }

    apex_nudge <- if (!is.null(offsets)) {
        0.05 * diff(range(dataset$intensity_plot, na.rm = TRUE)) /
            max(1L, length(offsets))
    } else {
        NULL
    }

    extra_layers <- list(
        highlight_peaks(
            dataset, supporting_datasets$detected_peaks, options,
            y_transform = transform_fun,
            offsets = offsets,
            stack_col = stack_col),
        highlight_apices(dataset, options, grouping_vars, nudge = apex_nudge),
        highlight_spectra_scans(datasets$spectra, options),
        rt_lines(options),
        legend_title(options),
        faceting(options, single),
        grid_layout(options, single)
    )
    extra_layers <- remove_null_elements(extra_layers)

    x_label <- ifelse(
        options$chromatograms$rt_unit == "minute",
        "RT (minutes)",
        "RT (seconds)")
    y_label <- ifelse(
        options$chromatograms$intensity_unit == "relative",
        "Relative intensity",
        "Intensity")

    # A stacked axis has no single meaning, so xcms drops it. Do the same.
    if (!is.null(offsets)) {
        y_label <- NULL
        extra_layers <- c(extra_layers, list(scale_y_continuous(breaks = NULL)))
    }

    p <- ggplot(
        data = dataset,
        mapping = build_aes(
            x = "rt_plot",
            y = "intensity_plot",
            options = options,
            group = "series_id")
    ) +
        geom_line(linetype = options$chromatograms$line_type) +
        labs(x = x_label, y = y_label) +
        scale_fill_discrete(guide = "none") + # Removes peak highlight legend
        scale_x_continuous(breaks = scales::pretty_breaks()) +
        theme_minimal() +
        extra_layers

    if (isTRUE(options$purity_overlay$show)) {
        p <- purity_overlay_layer(p, supporting_datasets$purity_scores, options)
    }

    # xcms returns the offsets so callers can draw their own reference lines.
    if (!is.null(offsets)) {
        attr(p, "stacked_offsets") <- offsets
    }

    return(p)
}
