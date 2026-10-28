.pdf_page_size <- function(page_size, orientation) {
    if (is.numeric(page_size)) {
        if (length(page_size) != 2) {
            stop(
                "lp_report: a numeric 'page_size' must be the width and ",
                "height in inches.", call. = FALSE)
        }
        return(page_size)
    }
    page_size <- match.arg(page_size, c("a4", "letter"))
    size <- switch(page_size, a4 = c(8.27, 11.69), letter = c(8.5, 11))
    if (orientation == "landscape") rev(size) else size
}

# The height of a row of a PDF table, in inches: a line of text and the
# padding of the table.
.pdf_row_height <- function(font_size) {
    font_size * 1.2 / 72 + 4 / 25.4
}

# The share of a PDF page the tables can take under the plots.
.pdf_table_share <- 0.45

.table_grob <- function(df, font_size, style) {
    cells <- function(colour, face, fill) {
        list(
            fg_params = list(
                hjust = 0, x = 0.05, col = colour, fontface = face),
            bg_params = list(fill = fill, col = NA))
    }
    grob <- gridExtra::tableGrob(
        df, rows = NULL,
        theme = gridExtra::ttheme_minimal(
            base_size = font_size,
            base_family = style$font_family,
            core = cells(style$text_colour, "plain", style$row_fills),
            colhead = cells(
                style$header_colour, style$header_face, style$header_fill)))
    if (style$rules) .add_rules(grob) else grob
}

# Rules above and below a table and under its header.
.add_rules <- function(grob) {
    rule <- function(y, lwd) {
        grid::segmentsGrob(
            x0 = 0, x1 = 1, y0 = y, y1 = y, gp = grid::gpar(lwd = lwd))
    }
    last <- nrow(grob)
    columns <- ncol(grob)
    grob <- gtable::gtable_add_grob(
        grob, rule(1, 1.5), t = 1, l = 1, r = columns, clip = "off")
    grob <- gtable::gtable_add_grob(
        grob, rule(0, 0.75), t = 1, l = 1, r = columns, clip = "off")
    gtable::gtable_add_grob(
        grob, rule(0, 1.5), t = last, l = 1, r = columns, clip = "off")
}

.pdf_text <- function(label, font_size, style, face = "plain",
                      colour = style$text_colour) {
    grid::textGrob(
        label, x = 0, hjust = 0,
        gp = grid::gpar(
            fontsize = font_size, fontface = face, col = colour,
            fontfamily = style$font_family))
}

# The grob at the top left of the space it is given, at its natural size.
.top_left <- function(grob) {
    if (inherits(grob, "gtable")) {
        grob$vp <- grid::viewport(
            x = 0, y = 1, just = c("left", "top"),
            width = sum(grob$widths), height = sum(grob$heights))
    } else {
        grob$y <- grid::unit(1, "npc") - grid::unit(0.5, "lines")
    }
    grob
}

# Draw a page of a PDF report within its margins.
.draw_pdf_page <- function(grob, margin = grid::unit(0.4, "in")) {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(
        width = grid::unit(1, "npc") - 2 * margin,
        height = grid::unit(1, "npc") - 2 * margin))
    grid::grid.draw(grob)
    grid::popViewport()
}

# The peaks table of a PDF page, cut to `max_rows` rows; the rest are
# summarised in a note.
.pdf_peaks_grob <- function(peaks, font_size, style, max_rows) {
    if (nrow(peaks) == 0) {
        return(.pdf_text("No detected peaks.", font_size, style, "italic"))
    }
    shown <- utils::head(peaks, max(1L, max_rows))
    grob <- .table_grob(shown, font_size, style)
    if (nrow(peaks) == nrow(shown)) {
        return(grob)
    }
    gridExtra::arrangeGrob(
        grob,
        .pdf_text(
            paste0(
                "... and ", nrow(peaks) - nrow(shown),
                " more peak(s); see the HTML report."),
            font_size, style, "italic"),
        ncol = 1,
        heights = grid::unit.c(sum(grob$heights), grid::unit(2, "lines")),
        widths = sum(grob$widths))
}

# The height of a table grob, or of a line of text.
.grob_height <- function(grob) {
    if (inherits(grob, "gtable")) sum(grob$heights) else grid::unit(2, "lines")
}

# The page of a PDF report: the heading, the plots, and the tables, laid out
# as one grob so that rendering errors surface before the file is written.
# The tables take their natural height, the plots the rest of the page.
.pdf_page_grob <- function(page, layout, size) {
    style <- layout$style
    font_size <- layout$font_size * 0.8
    single <- length(page$items) == 1
    # One feature's tables sit side by side; the wider tables of several
    # features are stacked, so the peaks get the rows the annotations leave.
    stacked <- !single

    budget <- floor(
        .pdf_table_share * (size[2] - 0.8) / .pdf_row_height(font_size))
    tables <- list()
    if (!is.null(page$annotations)) {
        tables <- c(tables, list(
            .table_grob(page$annotations, font_size, style)))
        if (stacked) {
            budget <- budget - nrow(page$annotations) - 1
        }
    }
    if (!is.null(page$peaks)) {
        tables <- c(tables, list(
            .pdf_peaks_grob(page$peaks, font_size, style, budget - 2)))
    }

    grobs <- list(patchwork::patchworkGrob(.page_figure(page, layout)))
    heights <- grid::unit(1, "null")
    if (length(tables) > 0) {
        table_heights <- do.call(grid::unit.c, lapply(tables, .grob_height))
        tables <- lapply(tables, .top_left)
        gap <- grid::unit(0.5, "lines")
        if (stacked || length(tables) == 1) {
            grobs <- c(grobs, list(gridExtra::arrangeGrob(
                grobs = tables, ncol = 1, heights = table_heights)))
            heights <- grid::unit.c(heights, sum(table_heights) + gap)
        } else {
            grobs <- c(grobs, list(gridExtra::arrangeGrob(
                grobs = tables, ncol = 2, widths = c(1, 2))))
            heights <- grid::unit.c(heights, max(table_heights) + gap)
        }
    }

    gridExtra::arrangeGrob(
        grobs = grobs,
        ncol = 1,
        heights = heights,
        top = .pdf_text(
            .page_heading(page), font_size * 1.5, style,
            style$heading_face, style$heading_colour))
}

.write_report_pdf <- function(rendered, file, title, layout, size) {
    index <- rendered$index
    style <- layout$style
    font_size <- layout$font_size * 0.8

    # The rows that fit within the margins, under the title and the header.
    rows_per_page <- max(
        5L, floor((size[2] - 1.6) / .pdf_row_height(font_size)) - 1L)
    chunks <- split(
        seq_len(nrow(index)), ceiling(seq_len(nrow(index)) / rows_per_page))
    if (length(chunks) == 0) {
        chunks <- list(integer())
    }

    shown <- .report_index_table(index, layout$by, layout$rt_unit)
    shown$page <- ifelse(
        is.na(index$page), "", as.character(length(chunks) + index$page))

    grDevices::pdf(file, width = size[1], height = size[2])
    on.exit(grDevices::dev.off())

    for (k in seq_along(chunks)) {
        heading <- if (k == 1) title else paste0(title, " (continued)")
        .draw_pdf_page(gridExtra::arrangeGrob(
            .top_left(if (length(chunks[[k]]) > 0) {
                .table_grob(
                    shown[chunks[[k]], , drop = FALSE], font_size, style)
            } else {
                .pdf_text("No features.", font_size, style, "italic")
            }),
            top = .pdf_text(
                heading, font_size * 1.8, style, style$title_face,
                style$heading_colour)))
    }

    for (page in rendered$pages) {
        .draw_pdf_page(page$figure)
    }
}
