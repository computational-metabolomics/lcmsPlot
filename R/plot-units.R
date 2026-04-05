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

highlight_peaks <- function(dataset, detected_peaks, options) {
    if (!options$chromatograms$highlight_peaks || nrow(detected_peaks) == 0) {
        return(NULL)
    }

    mode <- options$chromatograms$highlight_peaks_mode
    highlight_peaks_color <- options$chromatograms$highlight_peaks_color
    highlight_peaks_factor <- options$chromatograms$highlight_peaks_factor

    if (mode == "polygon") {
        plot_data <- build_ribbon_data(dataset, detected_peaks)
    } else {
        rt_scale <- if (options$chromatograms$rt_unit == "minute") 60 else 1
        plot_data <- detected_peaks |>
            mutate(
                rt_plot = .data$rt / rt_scale,
                rtmin_plot = .data$rtmin / rt_scale,
                rtmax_plot = .data$rtmax / rt_scale
            )
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
                       ymin = 0,
                       ymax = .data$intensity,
                       group = .data$peak_id,
                       fill = .data[[factor]],
                       colour = .data[[factor]]
                   ))))
               } else {
                   do.call(geom_ribbon, c(common, list(
                       mapping = aes(
                           ymin = 0,
                           ymax = .data$intensity,
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
                       ymin = 0,
                       ymax = .data$maxo,
                       fill = .data[[factor]],
                       colour = .data[[factor]]
                   ))))
               } else {
                   do.call(geom_rect, c(common, list(
                       mapping = aes(
                           xmin = .data$rtmin_plot,
                           xmax = .data$rtmax_plot,
                           ymin = 0,
                           ymax = .data$maxo
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
                       y = .data$maxo,
                       colour = .data[[factor]]
                   ))))
               } else {
                   do.call(geom_point, c(common, list(
                       mapping = aes(
                           x = .data$rt_plot,
                           y = .data$maxo
                       ),
                       colour = color
                   )))
               }
           },
           stop("Unknown highlight_peaks_mode: '", mode, "'. ",
                "Expected 'polygon', 'rectangle', or 'point'.")
    )
}

highlight_apices <- function(dataset, options, grouping_vars) {
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
        geom_text(
            data = top_peaks,
            aes(label = round(.data$rt_plot, 2)),
            nudge_y = 0.05 * max(dataset$intensity_plot),
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
