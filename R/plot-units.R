build_aes <- function(x, y, options, ...) {
    args <- list(...)
    p_aes <- list(x = sym(x), y = sym(y))

    if (!is.null(options$arrangement$group_by)) {
        p_aes$color <- sym(options$arrangement$group_by)
    }

    for (arg_name in names(args)) {
        p_aes[[arg_name]] <- sym(args[[arg_name]])
    }

    return(do.call(aes, p_aes))
}

highlight_peaks_aes <- function(ymax, options) {
    p_aes <- list(ymax = sym(ymax), ymin = 0)

    if (!is.null(options$arrangement$group_by)) {
        p_aes$fill <- sym(options$arrangement$group_by)
    }

    return(do.call(aes, p_aes))
}

#' Offsets that stack the series of a chromatogram panel up the y axis
#'
#' Reproduces the offsets computed by `xcms::plotChromatogramsOverlay()`: the
#' series are ordered by m/z and mapped linearly onto `stacked * y_limits`, so
#' the largest m/z sits at the top of the band.
#'
#' @param dataset A `data.frame` of chromatogram rows carrying `stack_col`, and
#' `feature_mz` when the input provides per-series m/z values.
#' @param stack_col A `character` naming the column that identifies a series.
#' @param stacked A `numeric` fraction of the intensity range to spread the
#' series over.
#' @param y_limits A length-2 `numeric` giving the transformed intensity range.
#' @return A named `numeric` vector of offsets, one per series key.
#' @keywords internal
stacked_offsets <- function(dataset, stack_col, stacked, y_limits) {
    keys <- unique(as.character(dataset[[stack_col]]))

    mzs <- if ("feature_mz" %in% names(dataset)) {
        vapply(keys, function(k) {
            mean(
                dataset$feature_mz[as.character(dataset[[stack_col]]) == k],
                na.rm = TRUE)
        }, numeric(1))
    } else {
        rep(NA_real_, length(keys))
    }

    # xcms falls back to the series position when no m/z is available. Series
    # that share an m/z need the same treatment, or they would all land on one
    # offset and nothing would appear stacked at all.
    if (anyNA(mzs) || any(!is.finite(mzs)) ||
        length(unique(mzs)) < length(keys)) {
        mzs <- as.numeric(rank(keys, ties.method = "first"))
    }

    band <- stacked * y_limits
    mz_range <- range(mzs)

    offsets <- if (diff(mz_range) == 0) {
        # A single series: nothing to spread out.
        rep(band[2L], length(keys))
    } else {
        band[2L] + (mzs - mz_range[2L]) * (band[2L] - band[1L]) / diff(mz_range)
    }

    stats::setNames(as.numeric(offsets), keys)
}

#' Look up the stacking offset for every row of a data frame
#'
#' @param df A `data.frame` carrying `stack_col`.
#' @param offsets A named `numeric` vector as returned by [stacked_offsets()],
#' or `NULL` when the plot is not stacked.
#' @param stack_col A `character` naming the series column, or `NULL`.
#' @return A `numeric` vector of offsets, zero where no offset applies.
#' @keywords internal
series_offset <- function(df, offsets, stack_col) {
    if (is.null(offsets) || is.null(stack_col) || !stack_col %in% names(df)) {
        return(rep(0, nrow(df)))
    }

    out <- unname(offsets[as.character(df[[stack_col]])])
    out[is.na(out)] <- 0
    out
}

highlight_peaks <- function(
    dataset,
    detected_peaks,
    options,
    y_transform = identity,
    offsets = NULL,
    stack_col = NULL
) {
    if (!options$chromatograms$highlight_peaks || nrow(detected_peaks) == 0) {
        return(NULL)
    }

    mode <- options$chromatograms$highlight_peaks_mode
    highlight_peaks_color <- options$chromatograms$highlight_peaks_color
    highlight_peaks_factor <- options$chromatograms$highlight_peaks_factor

    if (mode == "polygon") {
        # Rows are sliced out of `dataset`, whose intensity_plot already carries
        # the transform and the offset.
        plot_data <- build_ribbon_data(dataset, detected_peaks)
        plot_data$y_base <- series_offset(plot_data, offsets, stack_col)
    } else {
        rt_scale <- if (options$chromatograms$rt_unit == "minute") 60 else 1
        plot_data <- detected_peaks |>
            mutate(
                rt_plot = .data$rt / rt_scale,
                rtmin_plot = .data$rtmin / rt_scale,
                rtmax_plot = .data$rtmax / rt_scale
            )
        # `detected_peaks` is a different frame from `dataset`, so the apex
        # heights have to be transformed and offset here or the boxes detach
        # from the traces.
        plot_data$y_base <- series_offset(plot_data, offsets, stack_col)
        plot_data$maxo_plot <- y_transform(plot_data$maxo) + plot_data$y_base
    }

    build_geom(mode, plot_data, highlight_peaks_color, highlight_peaks_factor)
}

build_ribbon_data <- function(dataset, detected_peaks) {
    detected_peaks |>
        mutate(peak_id = row_number()) |>
        rowwise() |>
        do({
            peak <- .
            dataset |>
                filter(.data$sample_id == peak$sample_id) |>
                filter(.data$rt >= peak$rtmin, .data$rt <= peak$rtmax) |>
                mutate(peak_id = peak$peak_id)
        }) |>
        bind_rows()
}

build_geom <- function(mode, plot_data, color, factor) {
    has_color <- !is.null(color)

    switch(mode,
           polygon = {
               common <- list(data = plot_data, alpha = 0.3, linetype = 1)
               if (!has_color) {
                   do.call(geom_ribbon, c(common, list(mapping = aes(
                       ymin = .data$y_base,
                       ymax = .data$intensity_plot,
                       group = .data$peak_id,
                       fill = .data[[factor]],
                       colour = .data[[factor]]
                   ))))
               } else {
                   do.call(geom_ribbon, c(common, list(
                       mapping = aes(
                           ymin = .data$y_base,
                           ymax = .data$intensity_plot,
                           group = .data$peak_id
                       ),
                       fill = color,
                       colour = color
                   )))
               }
           },
           rectangle = {
               common <- list(data = plot_data, alpha = 0.2, inherit.aes = FALSE)
               if (!has_color) {
                   do.call(geom_rect, c(common, list(mapping = aes(
                       xmin = .data$rtmin_plot,
                       xmax = .data$rtmax_plot,
                       ymin = .data$y_base,
                       ymax = .data$maxo_plot,
                       fill = .data[[factor]],
                       colour = .data[[factor]]
                   ))))
               } else {
                   do.call(geom_rect, c(common, list(
                       mapping = aes(
                           xmin = .data$rtmin_plot,
                           xmax = .data$rtmax_plot,
                           ymin = .data$y_base,
                           ymax = .data$maxo_plot
                       ),
                       fill = color,
                       colour = color
                   )))
               }
           },
           point = {
               common <- list(data = plot_data, size = 2, inherit.aes = FALSE)
               if (!has_color) {
                   do.call(geom_point, c(common, list(mapping = aes(
                       x = .data$rt_plot,
                       y = .data$maxo_plot,
                       colour = .data[[factor]]
                   ))))
               } else {
                   do.call(geom_point, c(common, list(
                       mapping = aes(
                           x = .data$rt_plot,
                           y = .data$maxo_plot
                       ),
                       colour = color
                   )))
               }
           },
           stop("Unknown highlight_peaks_mode: '", mode, "'. ",
                "Expected 'polygon', 'rectangle', or 'point'.")
    )
}

highlight_apices <- function(dataset, options, grouping_vars, nudge = NULL) {
    highlight_apices_opts <- options$chromatograms$highlight_apices

    top_peaks <- dataset |>
        group_by(across(all_of(grouping_vars))) |>
        group_modify(~ {
            if (!is.null(highlight_apices_opts$column)) {
                hp <- unique(.x[[highlight_apices_opts$column]])
            } else {
                hp <- NA
            }

            if (!is.na(hp)) {
                # filter within ±5 rt and select max
                .x |>
                    filter(rt >= hp - 5, rt <= hp + 5) |>
                    slice_max(intensity, n = 1)
            } else if (!is.null(highlight_apices_opts$top_n)) {
                # take global maximum
                .x |> slice_max(intensity, n = highlight_apices_opts$top_n)
            } else {
                tibble()
            }
        })

    if (nrow(top_peaks) > 0) {
        # On a stacked plot the global maximum spans every series, so the caller
        # supplies a nudge scaled to one series' band instead.
        if (is.null(nudge)) {
            nudge <- 0.05 * max(dataset$intensity_plot)
        }

        geom_text(
            data = top_peaks,
            aes(label = round(.data$rt_plot, 2)),
            nudge_y = nudge,
            size = 3,
            color = "red"
        )
    } else {
        NULL
    }
}

highlight_spectra_scans <- function(dataset, options) {
    if (options$spectra$show) {
        geom_vline(
            data = dataset,
            aes(xintercept = .data$rt),
            color = "black",
            linetype = "dashed")
    } else {
        NULL
    }
}

rt_lines <- function(options) {
    lines <- lapply(options$rt_lines, function (rt_line_obj) {
        geom_vline(
            xintercept = rt_line_obj$intercept,
            color = rt_line_obj$color,
            linetype = rt_line_obj$line_type)
    })
}

legend_title <- function(options) {
    if (!is.null(options$labels$legend)) {
        return(guides(
            color = guide_legend(title = options$labels$legend),
            fill = guide_legend(title = options$labels$legend)
        ))
    } else {
        return(NULL)
    }
}

faceting <- function(options, single) {
    if (single & !is.null(options$facets$facets)) {
        if (options$facets$free_x & options$facets$free_y) {
            scales <- "free"
        } else if (options$facets$free_x) {
            scales <- "free_x"
        } else if (options$facets$free_y) {
            scales <- "free_y"
        } else {
            scales <- "fixed"
        }

        return(facet_wrap(
            as.formula(paste(
                "~",
                paste(options$facets$facets, collapse = "+"))),
            ncol = options$facets$ncol,
            nrow = options$facets$nrow,
            scales = scales))
    } else {
        return(NULL)
    }
}

grid_layout <- function(options, single) {
    rows <- options$grid$rows
    cols <- options$grid$cols

    if (single & (!is.null(rows) | !is.null(cols))) {
        facet_formula <- if (!is.null(rows) & !is.null(cols)) {
            as.formula(paste(rows, "~", cols))
        } else if (!is.null(rows)) {
            as.formula(paste(rows, "~ ."))
        } else if (!is.null(cols)) {
            as.formula(paste(". ~", cols))
        } else {
            NULL  # No faceting
        }

        if (options$grid$free_x && options$grid$free_y) {
            scales <- "free"
        } else if (options$grid$free_x) {
            scales <- "free_x"
        } else if (options$grid$free_y) {
            scales <- "free_y"
        } else {
            scales <- "fixed"
        }

        return(facet_grid(facet_formula, scales = scales))
    } else {
        return(NULL)
    }
}

title_tile <- function(text) {
    text_grob <- grid::textGrob(
        text,
        gp = grid::gpar(col = "black", fontsize = 9))
    return(text_grob)
}
