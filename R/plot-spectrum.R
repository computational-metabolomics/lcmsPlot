#' Plot spectra from one or more datasets
#'
#' `plot_spectrum()` generates a ggplot2 MS spectra plot.
#'
#' @param datasets A named `list` of data frames containing the primary datasets
#' to plot. For a spectra plot the used key is `spectra`.
#' @param supporting_datasets A `list` of supporting data frames, such as
#' detected peaks, used for highlighting features.
#' @param options A list of plot options, controlling units, faceting,
#' highlighting, and other visual parameters.
#' @param single A `logical` value that indicates whether it should treat
#' the plot as a single dataset variant, which can affect faceting
#' and layout behavior. Default is `FALSE`.
#' @return A `ggplot` object representing the RT difference plot.
#' @keywords internal
plot_spectrum <- function(
    datasets,
    supporting_datasets,
    options,
    single = FALSE
) {
    dataset <- datasets$spectra

    extra_layers <- list(
        legend_title(options),
        faceting(options, single)
    )
    extra_layers <- remove_null_elements(extra_layers)

    dataset_for_plot <- dataset |>
        mutate(intensity = case_when(
            reference == TRUE ~ -.data$intensity,
            TRUE ~ .data$intensity
        )) |>
        mutate(
            sample_id_rt = paste0(
                .data$sample_id, " RT: ", round(.data$rt, 3), " sec.",
                .precursor_label(dataset[["precursor_mz"]])
            )
        )

    top_peaks <- dataset_for_plot |>
        group_by(.data$sample_id_rt) |>
        slice_max(order_by = .data$intensity, n = 3)

    if (length(unique(dataset_for_plot$reference)) == 2) {
        p_aes <- build_aes(
            x = "mz",
            y = "intensity",
            options = options,
            color = "reference")
        colour_scales <- .palette_scale(options, "colour")
    } else {
        p_aes <- build_aes(
            x = "mz",
            y = "intensity",
            options = options)
        colour_scales <- .series_scales(options)
    }

    # Determine dynamic y-axis limits
    y_max <- ceiling(max(dataset_for_plot$intensity) / 10) * 10
    y_min <- floor(min(dataset_for_plot$intensity) / 10) * 10

    p <- ggplot(
        data = dataset_for_plot,
        mapping = p_aes
    ) +
        geom_segment(
            aes(xend = .data$mz, yend = 0),
            linewidth = .line_width(options),
            show.legend = single) +
        geom_text(
            data = top_peaks,
            aes(label = round(.data$mz, 4)),
            nudge_y = 0.05 * y_max,
            size = options$spectra$peak_label_size,
            color = "red"
        ) +
        labs(x = "m/z", y = "Relative intensity (%)") +
        scale_x_continuous(
            breaks = scales::pretty_breaks(n = options$spectra$mz_breaks_n)) +
        scale_y_continuous(
            limits = c(y_min * 1.1, y_max * 1.1),
            breaks = seq(y_min, y_max, by = options$spectra$intensity_breaks_by),
            labels = function(x) abs(x)
        ) +
        colour_scales +
        .base_theme(options) +
        extra_layers

    opt <- options$spectra
    precursors <- .precursor_markers(
        dataset_for_plot, opt$precursor_marker_position, opt$ppm)
    if (nrow(precursors) > 0) {
        p <- p + geom_point(
            data = precursors,
            aes(x = .data$precursor_mz, y = .data$marker_y),
            inherit.aes = FALSE,
            shape = opt$precursor_marker_shape,
            size = 2.5,
            color = opt$precursor_marker_color,
            fill = opt$precursor_marker_color)
    }

    if (is.null(options$facets$facets) && options$spectra$auto_facet) {
        p <- p + facet_wrap(~ sample_id_rt, ncol = 1)
    }

    return(p)
}

# The precursor m/z appended to the facet label of an MS2 spectrum.
.precursor_label <- function(precursor_mz) {
    if (is.null(precursor_mz)) {
        return("")
    }
    ifelse(
        is.na(precursor_mz),
        "",
        paste0(" Precursor m/z: ", round(precursor_mz, 4)))
}

# One row per spectrum with a precursor, carrying the faceting columns and the
# height of its marker: on the m/z axis (`position = "bottom"`), or on top of
# the precursor's peak (`"top"`), the highest within `ppm` of the precursor
# m/z, falling back to the m/z axis when the spectrum has no such peak.
.precursor_markers <- function(dataset, position, ppm) {
    if (!"precursor_mz" %in% colnames(dataset)) {
        return(tibble())
    }
    dataset |>
        filter(!.data$reference, !is.na(.data$precursor_mz)) |>
        group_by(across(all_of(c("sample_id_rt", "precursor_mz")))) |>
        mutate(marker_y = switch(
            position,
            bottom = 0,
            top = .precursor_peak_height(
                .data$mz, .data$intensity, .data$precursor_mz[1], ppm))) |>
        slice(1) |>
        ungroup()
}

.precursor_peak_height <- function(mz, intensity, precursor_mz, ppm) {
    near <- abs(mz - precursor_mz) <= precursor_mz * ppm / 1e6
    if (any(near)) max(intensity[near]) else 0
}
