# The heading of a feature: the value of `by`, with the name of the feature
# when it is another column.
.report_heading <- function(row, by) {
    key <- row[[by]][1]
    name <- if ("name" %in% colnames(row)) row$name[1] else NA
    if (by == "name" || is.na(name) || identical(as.character(name), key)) {
        return(key)
    }
    paste0(key, ": ", name)
}

# The heading of a page: its feature, or the features it gathers.
.page_heading <- function(page) {
    if (length(page$items) == 1) {
        return(page$items[[1]]$heading)
    }
    keys <- vapply(page$items, `[[`, character(1), "key")
    label <- paste(keys, collapse = ", ")
    if (nchar(label) > 80) {
        label <- paste0(keys[1], ", ..., ", keys[length(keys)])
    }
    label
}

# The detected peak of a feature in every sample.
.report_peaks <- function(obj) {
    peaks <- obj@data@detected_peaks
    if (nrow(peaks) == 0) {
        return(tibble())
    }
    peaks <- .join_sample_data(peaks, .sample_table(obj@data)) |>
        arrange(.data$sample_index)
    several <- "feature_id" %in% colnames(peaks) &&
        length(unique(peaks$feature_id)) > 1
    columns <- c(
        "sample_id", if (several) "feature_id",
        "mz", "rt", "rtmin", "rtmax", "into", "maxo")
    peaks |> select(any_of(columns))
}

# The tables of a page, formatted for display, or NULL when their section is
# off. A page of one feature lists its annotations as fields and values; a
# page of several, one feature per row, and its peaks carry their feature.
.page_tables <- function(items, layout) {
    by <- layout$by
    single <- length(items) == 1

    annotations <- NULL
    if ("annotations" %in% layout$sections) {
        rows <- bind_rows(lapply(items, `[[`, "annotations"))
        annotations <- .report_index_table(rows, by, layout$rt_unit)
        annotations$status <- NULL
        annotations$page <- NULL
        if (single) {
            annotations <- tibble(
                field = colnames(annotations),
                value = unlist(annotations, use.names = FALSE))
        }
    }

    peaks <- NULL
    if ("peaks" %in% layout$sections) {
        peaks <- bind_rows(lapply(items, function(item) {
            if (single || nrow(item$peaks) == 0) {
                return(item$peaks)
            }
            bind_cols(tibble(!!by := item$key), item$peaks)
        }))
        peaks <- .format_table(peaks, layout$rt_unit)
    }

    list(annotations = annotations, peaks = peaks)
}

# The table with its numbers formatted for display: m/z to 4 decimals, rt to
# 2 in the unit of the plots (named in the rt column headers), and other
# numbers to 4 significant digits. Missing values are blank.
.format_table <- function(df, rt_unit = "second") {
    df <- as.data.frame(df)
    rt_columns <- intersect(c("rt", "rtmin", "rtmax"), colnames(df))
    for (column in colnames(df)) {
        x <- df[[column]]
        formatted <- if (column %in% c("mz", "mzmin", "mzmax") &&
                         is.numeric(x)) {
            sprintf("%.4f", x)
        } else if (column %in% rt_columns && is.numeric(x)) {
            sprintf("%.2f", if (rt_unit == "minute") x / 60 else x)
        } else if (is.numeric(x) && !is.integer(x)) {
            formatC(x, digits = 4, format = "g")
        } else {
            as.character(x)
        }
        formatted[is.na(x)] <- ""
        df[[column]] <- formatted
    }
    suffix <- if (rt_unit == "minute") " (min)" else " (s)"
    colnames(df)[match(rt_columns, colnames(df))] <- paste0(rt_columns, suffix)
    df
}

.report_rt_unit <- function(options) {
    options$chromatograms$rt_unit %||% "second"
}

# The index columns shown in a report, with the feature IDs listed when a
# value of `by` can cover several.
.report_index_table <- function(index, by, rt_unit) {
    shown <- index
    shown$feature_ids <- if (by == "feature_id") {
        NULL
    } else {
        vapply(shown$feature_ids, paste, character(1), collapse = "; ")
    }
    shown$page <- NULL
    .format_table(shown, rt_unit)
}

# The figure of a page: the plots of its features, each above its MS2
# spectra when there are any, in a grid of `layout$ncol` columns. With
# several features, each plot is titled with its feature.
.page_figure <- function(page, layout) {
    style <- layout$style
    single <- length(page$items) == 1
    panels <- lapply(page$items, function(item) {
        parts <- list(patchwork::wrap_elements(
            full = patchwork::patchworkGrob(item$plot)))
        if (!is.null(item$ms2_plot)) {
            parts <- c(parts, list(patchwork::wrap_elements(
                full = patchwork::patchworkGrob(item$ms2_plot))))
        }
        panel <- patchwork::wrap_plots(
            parts, ncol = 1, heights = c(3, 2)[seq_along(parts)])
        if (single) {
            return(panel)
        }
        patchwork::wrap_elements(full = patchwork::patchworkGrob(panel)) +
            labs(title = item$heading) +
            theme(plot.title = element_text(
                family = style$font_family,
                face = style$heading_face,
                colour = style$heading_colour))
    })
    if (single) {
        return(panels[[1]])
    }
    patchwork::wrap_plots(panels, ncol = min(layout$ncol, length(panels)))
}

# The features of a page with no MS2 spectra, when the section is on.
.missing_ms2 <- function(page) {
    has_ms2 <- any(vapply(page$items, `[[`, logical(1), "has_ms2"))
    if (!has_ms2) {
        return(character())
    }
    missing <- !vapply(page$items, `[[`, logical(1), "ms2_found")
    vapply(page$items[missing], `[[`, character(1), "key")
}
