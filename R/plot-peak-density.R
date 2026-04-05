#' Plot a peak density plot
#'
#' `plot_peak_density()` generates a ggplot2 peak density plot that mirrors
#' `xcms::plotChromPeakDensity()`. For each m/z feature:
#' - The x axis shows retention time.
#' - The y axis shows sample indices (1 to n), positioned within the density
#'   range so that sample rows are evenly spaced from 0 to `max(density)`.
#' - Detected peaks are plotted as coloured points at each sample's y position.
#' - The kernel density estimate of peak apex RTs is drawn as a line.
#' - When feature group simulation is enabled, semi-transparent rectangles
#'   mark feature groups that pass the `min_fraction` / `min_samples`
#'   threshold.
#'
#' @param datasets A named `list` of data frames. The used key is
#' `peak_density`.
#' @param supporting_datasets A `list` of supporting data frames. The used key
#' is `detected_peaks`, which provides per-sample peak positions.
#' @param options A list of plot options.
#' @param single A `logical` value indicating whether to treat the plot as a
#' single dataset variant. Default is `FALSE`.
#' @return A `ggplot` object representing the peak density plot.
#' @keywords internal
plot_peak_density <- function(
    datasets,
    supporting_datasets,
    options,
    single = FALSE
) {
    peak_density_data <- datasets$peak_density
    detected_peaks <- supporting_datasets$detected_peaks
    features <- options$peak_density$features

    # Split density curve rows from feature-group rectangle rows
    density_data <- peak_density_data |>
        filter(.data$data_type == "density")

    rect_data <- peak_density_data |>
        filter(.data$data_type == "rect")

    # n_samples: total number of samples (maximum sample index in the study)
    n_samples <- if (nrow(detected_peaks) > 0) {
        max(detected_peaks$sample_index, na.rm = TRUE)
    } else {
        0L
    }

    # ypos: evenly space sample indices within [0, global_max_density].
    # Peak for sample k is drawn at ypos[k], matching xcms.
    global_max_density <- max(density_data$density, na.rm = TRUE)
    ypos <- seq(from = 0, to = global_max_density, length.out = n_samples)

    rt_unit <- options$peak_density$rt_unit
    rt_scale <- if (rt_unit == "minute") 60 else 1
    x_label <- if (rt_unit == "minute") "RT (minutes)" else "RT (seconds)"

    density_data <- density_data |>
        mutate(rt_plot = .data$rt / rt_scale)

    # Feature label for auto-faceting: "m/z [mzmin - mzmax]"
    make_label <- function(mzmin, mzmax) {
        paste0("m/z [", round(mzmin, 4), " \u2013 ", round(mzmax, 4), "]")
    }

    density_data <- density_data |>
        mutate(feature_label = make_label(.data$mzmin, .data$mzmax))

    fid_to_label <- density_data |>
        distinct(.data$feature_metadata_id, .data$feature_label)

    # Per-feature max density (for rect ymax when faceting with free scales)
    feat_max_density <- density_data |>
        group_by(.data$feature_metadata_id) |>
        summarize(feat_max_density = max(.data$density), .groups = "drop")

    # Assign y positions and feature labels to detected peaks per feature
    peaks_by_feature <- lapply(seq_len(nrow(features)), function(i) {
        feat <- features[i, ]

        pks <- detected_peaks |>
            filter(
                .data$mz >= feat[["mzmin"]],
                .data$mz <= feat[["mzmax"]]
            )

        has_rt <- !is.na(feat[["rtmin"]]) && !is.na(feat[["rtmax"]])
        if (has_rt) {
            pks <- pks |>
                filter(
                    .data$rt >= feat[["rtmin"]],
                    .data$rt <= feat[["rtmax"]]
                )
        }

        if (nrow(pks) == 0) return(NULL)

        pks |>
            mutate(
                rt_plot = .data$rt / rt_scale,
                y_pos = ypos[.data$sample_index],
                feature_metadata_id = as.integer(i)
            )
    })
    peaks_by_feature <- bind_rows(remove_null_elements(peaks_by_feature))

    if (nrow(peaks_by_feature) > 0) {
        peaks_by_feature <- peaks_by_feature |>
            left_join(fid_to_label, by = "feature_metadata_id")
    }

    if (nrow(rect_data) > 0) {
        rect_data <- rect_data |>
            left_join(fid_to_label, by = "feature_metadata_id") |>
            left_join(feat_max_density, by = "feature_metadata_id") |>
            mutate(
                rtmin_plot = .data$rtmin / rt_scale,
                rtmax_plot = .data$rtmax / rt_scale
            )
    }

    extra_layers <- list(
        legend_title(options),
        faceting(options, single),
        grid_layout(options, single)
    )
    extra_layers <- remove_null_elements(extra_layers)

    p <- ggplot(density_data, aes(x = .data$rt_plot, y = .data$density)) +
        geom_line(colour = "black") +
        scale_y_continuous(
            breaks = ypos,
            labels = seq_len(n_samples)
        ) +
        labs(x = x_label, y = "sample", colour = NULL) +
        theme_minimal() +
        extra_layers

    if (nrow(peaks_by_feature) > 0) {
        p <- p + geom_point(
            data = peaks_by_feature,
            mapping = aes(
                x = .data$rt_plot,
                y = .data$y_pos,
                colour = .data$sample_id
            ),
            inherit.aes = FALSE
        )
    }

    if (nrow(rect_data) > 0) {
        p <- p + geom_rect(
            data = rect_data,
            mapping = aes(
                xmin = .data$rtmin_plot,
                xmax = .data$rtmax_plot,
                ymin = 0,
                ymax = .data$feat_max_density
            ),
            fill = "#00000020",
            colour = "#00000040",
            inherit.aes = FALSE
        )
    }

    n_features <- length(unique(density_data$feature_metadata_id))
    auto_facet <- n_features > 1 &&
        single &&
        is.null(options$facets$facets) &&
        is.null(options$grid$rows) &&
        is.null(options$grid$cols)

    if (auto_facet) {
        p <- p + facet_wrap(~ feature_label, scales = "free")
    }

    return(p)
}
