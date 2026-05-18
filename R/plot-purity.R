#' Add purity score overlay to a chromatogram plot
#'
#' Adds a `geom_point` layer to a chromatogram ggplot where each point marks
#' the retention time of an MS/MS acquisition event, coloured by
#' interpolated precursor ion purity (`inPurity`).
#'
#' @param p A `ggplot` object (a chromatogram).
#' @param purity_scores A `data.frame` with columns `rt`, `in_purity`,
#'   `metadata_index`, and any additional sample columns merged in by the
#'   render pipeline.
#' @param options A `list` of plot options (reads `options$purity_overlay`).
#' @return The modified `ggplot` object.
#' @keywords internal
purity_overlay_layer <- function(p, purity_scores, options) {
    if (is.null(purity_scores) || nrow(purity_scores) == 0) {
        return(p)
    }

    rt_scale <- if (
        !is.null(options$chromatograms$rt_unit) &&
        options$chromatograms$rt_unit == "minute"
    ) 60 else 1

    overlay_data <- purity_scores |>
        dplyr::mutate(rt_plot = .data$rt / rt_scale)

    threshold <- options$purity_overlay$threshold

    midpoint <- if (!is.null(threshold)) threshold else 0.5

    purity_layers <- list(
        geom_point(
            data = overlay_data,
            aes(
                x = .data$rt_plot,
                y = -Inf,
                color = .data$in_purity
            ),
            size = options$purity_overlay$point_size,
            shape = 17,
            inherit.aes = FALSE
        ),
        scale_color_gradient2(
            name = "inPurity",
            low = "#d73027", mid = "#fee090", high = "#1a9850",
            midpoint = midpoint,
            limits = c(0, 1)
        )
    )

    if (requireNamespace("ggnewscale", quietly = TRUE)) {
        p + ggnewscale::new_scale_color() + purity_layers
    } else {
        p + purity_layers
    }
}

#' Plot precursor ion purity scores over retention time
#'
#' Generates a scatter plot of `inPurity` (y-axis) versus retention time
#' (x-axis) for each sample, with an optional horizontal threshold line.
#'
#' @param datasets A named `list` containing a `purity_timeline` element.
#' @param supporting_datasets Unused; present for API consistency.
#' @param options A `list` of plot options.
#' @param single Unused; present for API consistency.
#' @return A `ggplot` object.
#' @keywords internal
plot_purity_timeline <- function(
    datasets,
    supporting_datasets,
    options,
    single = FALSE
) {
    dataset <- datasets$purity_timeline

    rt_scale <- if (
        !is.null(options$chromatograms$rt_unit) &&
        options$chromatograms$rt_unit == "minute"
    ) 60 else 1

    dataset <- dataset |>
        dplyr::mutate(rt_plot = .data$rt / rt_scale)

    x_label <- if (rt_scale == 60) "RT (minutes)" else "RT (seconds)"
    threshold <- options$purity_timeline$threshold

    extra_layers <- list(
        legend_title(options),
        if (!is.null(threshold)) {
            geom_hline(
                yintercept = threshold,
                linetype = "dashed",
                color = "red",
                alpha = 0.7
            )
        }
    )
    extra_layers <- remove_null_elements(extra_layers)

    ggplot(
        data = dataset,
        mapping = build_aes(
            x = "rt_plot",
            y = "in_purity",
            options = options,
            group = "sample_id"
        )
    ) +
        geom_point(alpha = 0.6) +
        scale_color_gradient2(
            name = "inPurity",
            low = "#d73027", mid = "#fee090", high = "#1a9850",
            midpoint = if (!is.null(threshold)) threshold else 0.5,
            limits = c(0, 1)
        ) +
        labs(x = x_label, y = "Precursor ion purity") +
        scale_y_continuous(limits = c(0, 1)) +
        theme_minimal() +
        extra_layers
}

#' Plot the distribution of precursor ion purity scores per sample
#'
#' Generates a violin or boxplot of `inPurity` scores grouped by sample,
#' with an optional threshold line.
#'
#' @param datasets A named `list` containing a `purity_distribution` element.
#' @param supporting_datasets Unused; present for API consistency.
#' @param options A `list` of plot options.
#' @param single Unused; present for API consistency.
#' @return A `ggplot` object.
#' @keywords internal
plot_purity_distribution <- function(
    datasets,
    supporting_datasets,
    options,
    single = FALSE
) {
    dataset <- datasets$purity_distribution

    threshold <- options$purity_distribution$threshold
    geom_type <- options$purity_distribution$type
    if (is.null(geom_type)) geom_type <- "violin"
    geom_func <- get(paste0("geom_", geom_type), asNamespace("ggplot2"))

    extra_layers <- list(
        legend_title(options),
        if (!is.null(threshold)) {
            geom_hline(
                yintercept = threshold,
                linetype = "dashed",
                color = "red",
                alpha = 0.7
            )
        }
    )
    extra_layers <- remove_null_elements(extra_layers)

    ggplot(
        data = dataset,
        mapping = build_aes(
            x = "sample_id",
            y = "in_purity",
            options = options
        )
    ) +
        geom_func() +
        labs(x = "Sample", y = "Precursor ion purity") +
        scale_y_continuous(limits = c(0, 1)) +
        theme_minimal() +
        extra_layers
}

#' Plot an MS1 spectrum with isolation window annotation
#'
#' Extends the standard spectrum plot by adding a semi-transparent rectangle
#' spanning the isolation window and a dashed vertical line at the precursor
#' m/z. Requires `in_purity`, `precursor_mz`, and
#' `isolation_window_half_width` columns in the dataset (merged from
#' `feature_metadata` by the render pipeline).
#'
#' @param datasets A named `list` containing a `spectra` element.
#' @param supporting_datasets Unused; present for API consistency.
#' @param options A `list` of plot options.
#' @param single A `logical` value passed through to the base spectrum plot.
#' @return A `ggplot` object.
#' @keywords internal
plot_isolation_window <- function(
    datasets,
    supporting_datasets,
    options,
    single = FALSE
) {
    p <- plot_spectrum(datasets, supporting_datasets, options, single)

    dataset <- datasets$spectra

    has_window_cols <- all(
        c("precursor_mz", "isolation_window_half_width") %in% colnames(dataset)
    )

    if (!has_window_cols) {
        return(p)
    }

    window_data <- dataset |>
        dplyr::distinct(
            .data$feature_metadata_id,
            .keep_all = TRUE
        ) |>
        dplyr::select(
            "feature_metadata_id",
            "precursor_mz",
            "isolation_window_half_width",
            dplyr::any_of("in_purity")
        )

    p <- p +
        geom_rect(
            data = window_data,
            aes(
                xmin = .data$precursor_mz - .data$isolation_window_half_width,
                xmax = .data$precursor_mz + .data$isolation_window_half_width,
                ymin = -Inf,
                ymax = Inf
            ),
            alpha = 0.12,
            fill = "#4575b4",
            inherit.aes = FALSE
        ) +
        geom_vline(
            data = window_data,
            aes(xintercept = .data$precursor_mz),
            color = "#d73027",
            linetype = "dashed",
            alpha = 0.8,
            inherit.aes = FALSE
        )

    if ("in_purity" %in% colnames(window_data)) {
        purity_val <- round(mean(window_data$in_purity, na.rm = TRUE), 3)
        p <- p + labs(subtitle = paste0("inPurity = ", purity_val))
    }

    zoom_factor <- options$isolation_window$zoom_factor %||% 3
    half <- max(window_data$isolation_window_half_width)
    p <- p + coord_cartesian(xlim = c(
        min(window_data$precursor_mz) - zoom_factor * half,
        max(window_data$precursor_mz) + zoom_factor * half
    ))

    return(p)
}
