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
                .data$sample_id, " RT: ", round(.data$rt, 3), " sec."
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
    } else {
        p_aes <- build_aes(
            x = "mz",
            y = "intensity",
            options = options)
    }

    # Determine dynamic y-axis limits
    y_max <- max(dataset_for_plot$intensity)

    # Round limits to nearest 10 for clean breaks
    y_max <- ceiling(y_max / 10) * 10

    p <- ggplot(
        data = dataset_for_plot,
        mapping = p_aes
    ) +
        geom_segment(aes(xend = .data$mz, yend = 0), show.legend = single) +
        geom_text(
            data = top_peaks,
            aes(label = round(.data$mz, 4)),
            nudge_y = 0.05 * max(dataset_for_plot$intensity),
            size = 3,
            color = "red"
        ) +
        labs(x = "m/z", y = "Relative intensity (%)") +
        expand_limits(y = max(dataset_for_plot$intensity) * 1.1) +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 20)) +
        scale_y_continuous(
            limits = c(0, y_max * 1.1),
            breaks = seq(0, y_max, by = 20),
            labels = function(x) abs(x)
        ) +
        theme_minimal() +
        extra_layers

    if (is.null(options$facets$facets)) {
        p <- p + facet_wrap(~ sample_id_rt, ncol = 1)
    }

    return(p)
}
