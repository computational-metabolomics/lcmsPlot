#' Write a report with one page per feature
#'
#' `lp_report()` turns a plot pipeline into a report of the features, compound
#' names, or any other annotation of the chromatograms. Each feature is drawn
#' with the pipeline restricted to it (see [subset_features()]), with a table
#' of its annotations (`name`, `formula`, `adduct`, `mz`, `rt`, ...) and of its
#' detected peak in every sample (`rt`, `into`, `maxo`). By default each
#' feature has a page of its own; `features_per_page` tiles several on a page.
#'
#' The data is extracted once, when the pipeline is built: the report only
#' narrows it down for each feature, so a report over many features does not
#' read the raw files again.
#'
#' The HTML report is a single self-contained file: it opens with an index of
#' the features, which can be searched and sorted and links to their pages,
#' and embeds the plots as PNG images. The PDF report is one multi-page file,
#' with the index on the first pages and the feature pages after it.
#'
#' Values of `features` that match no data are listed in the index as
#' `"no data"`, and the features of a page that fails to render as `"error"`,
#' rather than stopping the report.
#'
#' By default the report records how it was made: a closing provenance
#' section (the last page of a PDF report) lists the versions of lcmsPlot, R
#' and Bioconductor, the layers of the pipeline with their arguments, and the
#' input files. The full record, with the version of every loaded package, is
#' written next to the report as JSON (`features.html` gives
#' `features.html.provenance.json`), and embedded in the HTML report. See
#' [lp_save()].
#'
#' @param object An instance of class `lcmsPlotClass`, with the chromatograms
#' of its features (see `lp_chromatogram(features = )` or the layers of the
#' compound sources, such as [lp_compound_discoverer()]).
#' @param file The path of the report. Its extension (`.html` or `.pdf`) sets
#' the format unless `format` is given.
#' @param by The annotation column that defines a feature of the report:
#' `"feature_id"` (the default), `"name"`, or any other column of the
#' chromatograms' `chromData()`. A feature of the report covers every feature
#' of the data that shares a value of `by`.
#' @param format The report format, `"html"` or `"pdf"`. Inferred from the
#' extension of `file` when `NULL`.
#' @param features The values of `by` to report, in the order of the pages.
#' When `NULL`, every value found in the data.
#' @param sections The sections drawn under the plots: `"annotations"` (the
#' table of the feature's annotations) and `"peaks"` (the table of its
#' detected peaks), both by default, and optionally `"ms2"`, the MS2 spectra
#' of the feature by precursor (see `lp_spectra(mode = "precursor")`),
#' extracted once for the whole report and drawn under its plot. `"ms2"` is
#' ignored when the pipeline already plots them.
#' @param features_per_page The number of features drawn on a page, in a grid
#' of `ncol` columns. With several features, each plot is titled with its
#' feature, and the tables list the features of the page.
#' @param ncol The number of columns of the grid of features on a page. By
#' default, up to 2.
#' @param theme The look of the report: `"plain"` (the default, unstyled),
#' `"clean"` (light cards, shaded table headers and status badges), or
#' `"journal"` (serif type and ruled tables, as in print). The plots keep the
#' theme of the pipeline (see [lp_theme()]).
#' @param title The title of the report.
#' @param width The width of the plots of a page of the HTML report, in inches.
#' @param height The height of a row of plots of the HTML report, in inches.
#' @param dpi The resolution of the plots of the HTML report.
#' @param page_size The page size of the PDF report: `"a4"`, `"letter"`, or a
#' `numeric` vector of the width and height in inches.
#' @param orientation The page orientation of the PDF report, `"landscape"`
#' or `"portrait"`. Ignored when `page_size` is `numeric`.
#' @param progress A `logical` value indicating whether to show a progress bar
#' while the pages are rendered.
#' @param provenance A `logical` value indicating whether to record the
#' provenance of the report.
#' @param checksums A `logical` value indicating whether to record the MD5
#' checksum of every input file in the provenance. Reading large raw files can
#' take a while, so it is off by default.
#' @return The index of the report (invisibly): a `tibble` with one row per
#' value of `by`, its annotations, the `feature_ids` it covers, the number of
#' samples with a detected peak (`n_samples`), its `status` (`"ok"`,
#' `"no data"`, or `"error"`), and the number of its `page` (`NA` when it has
#' none).
#' @export
#' @examples
#' raw_files <- dir(
#'    system.file("cdf", package = "faahKO"),
#'    full.names = TRUE,
#'    recursive = TRUE)[1:2]
#'
#' features <- data.frame(
#'     feature_id = c("A", "B"),
#'     mzmin = c(334.9, 344.0),
#'     mzmax = c(335.1, 344.2),
#'     rtmin = c(2700, 3400),
#'     rtmax = c(2900, 3500))
#'
#' p <- lcmsPlot(raw_files) +
#'     lp_chromatogram(features = features) +
#'     lp_arrange(group_by = "sample_id")
#'
#' if (requireNamespace("htmltools", quietly = TRUE)) {
#'     lp_report(p, tempfile(fileext = ".html"))
#'
#'     ## Both features on one page, without the peak tables
#'     lp_report(
#'         p, tempfile(fileext = ".html"),
#'         sections = "annotations",
#'         features_per_page = 2,
#'         theme = "clean")
#' }
lp_report <- function(
    object,
    file,
    by = "feature_id",
    format = NULL,
    features = NULL,
    sections = c("annotations", "peaks"),
    features_per_page = 1,
    ncol = NULL,
    theme = "plain",
    title = "lcmsPlot feature report",
    width = 10,
    height = 6,
    dpi = 150,
    page_size = "a4",
    orientation = "landscape",
    progress = interactive(),
    provenance = TRUE,
    checksums = FALSE
) {
    if (!is(object, "lcmsPlotClass") || !.has_data(object)) {
        stop(
            "lp_report: 'object' must be an lcmsPlot object with data; ",
            "add a layer such as lp_chromatogram(features = ) first.",
            call. = FALSE)
    }
    format <- format %||% tolower(tools::file_ext(file))
    if (!format %in% c("html", "pdf")) {
        stop(
            "lp_report: unknown format '", format, "'; use a '.html' or ",
            "'.pdf' file, or set 'format'.", call. = FALSE)
    }
    rlang::check_installed(
        if (format == "html") "htmltools" else c("gridExtra", "gtable"),
        reason = paste0("to write ", toupper(format), " reports."))
    if (length(sections) > 0) {
        sections <- match.arg(
            sections, c("annotations", "peaks", "ms2"), several.ok = TRUE)
    }
    if (!is.numeric(features_per_page) || length(features_per_page) != 1 ||
        features_per_page < 1) {
        stop(
            "lp_report: 'features_per_page' must be a number of at least 1.",
            call. = FALSE)
    }
    features_per_page <- as.integer(features_per_page)
    ncol <- as.integer(ncol %||% min(features_per_page, 2L))
    theme <- match.arg(theme, names(.report_styles))
    orientation <- match.arg(orientation, c("landscape", "portrait"))

    layout <- list(
        by = by,
        sections = sections,
        per_page = features_per_page,
        ncol = ncol,
        style = .report_styles[[theme]],
        font_size = .theme_sizes(object@options)$base_size,
        rt_unit = .report_rt_unit(object@options))

    index <- .report_index(object, by, features)
    ms2 <- if ("ms2" %in% sections) .report_ms2(object) else NULL

    render <- if (format == "html") {
        function(page) .html_page_figure(page, layout, width, height, dpi)
    } else {
        size <- .pdf_page_size(page_size, orientation)
        function(page) .pdf_page_grob(page, layout, size)
    }
    rendered <- .report_pages(object, index, layout, ms2, render, progress)

    prov <- if (isTRUE(provenance)) {
        .provenance(object, checksums, output = list(
            file = basename(file), format = format, by = by,
            features = nrow(rendered$index),
            pages = length(rendered$pages)))
    }
    if (format == "html") {
        .write_report_html(rendered, file, title, layout, prov)
    } else {
        .write_report_pdf(rendered, file, title, layout, size, prov)
    }
    if (!is.null(prov)) {
        .write_provenance(prov, file)
    }

    invisible(rendered$index)
}

#' The MS2 spectra of every feature, extracted once for a report
#'
#' @param object An instance of class `lcmsPlotClass`.
#' @return A copy of `object` that plots only the MS2 spectra of its features
#' by precursor, or `NULL` when the pipeline already plots them or there are
#' none.
#' @keywords internal
.report_ms2 <- function(object) {
    spectra_opts <- object@options$spectra
    if (isTRUE(spectra_opts$show) &&
        identical(spectra_opts$mode, "precursor")) {
        message(
            "lp_report: the pipeline already plots the MS2 spectra by ",
            "precursor; the 'ms2' section is skipped.")
        return(NULL)
    }

    ms2 <- tryCatch(
        lp_spectra(mode = "precursor")(object, FALSE),
        error = function(e) {
            warning(
                "lp_report: the 'ms2' section is skipped: ",
                conditionMessage(e), call. = FALSE)
            NULL
        })
    if (is.null(ms2)) {
        return(NULL)
    }
    if (length(ms2@data@spectra) == 0) {
        warning(
            "lp_report: no MS2 spectra match the features; the 'ms2' ",
            "section is skipped.", call. = FALSE)
        return(NULL)
    }

    for (name in setdiff(names(ms2@options), "spectra")) {
        if (is.list(ms2@options[[name]]) &&
            !is.null(ms2@options[[name]]$show)) {
            ms2@options[[name]]$show <- FALSE
        }
    }
    ms2
}

#' Render the pages of a report
#'
#' Narrows the data down to each row of the index with [subset_features()],
#' by feature ID so that every dataset carrying one is narrowed, gathers the
#' features into pages of `layout$per_page`, and renders each page with
#' `render` as soon as it is full. A feature with no data is marked
#' `"no data"` in the index, and the features of a page that fails `"error"`.
#'
#' @param object An instance of class `lcmsPlotClass`.
#' @param index The index of the report, from [.report_index()].
#' @param layout A `list` of the layout settings of the report.
#' @param ms2 The MS2 object from [.report_ms2()], or `NULL`.
#' @param render A function of a page that returns its rendered figure.
#' @param progress A `logical` value indicating whether to show a progress
#' bar.
#' @return A `list` with the `index`, updated with the `status` and `page` of
#' every row, and the rendered `pages`.
#' @keywords internal
.report_pages <- function(object, index, layout, ms2, render, progress) {
    by <- layout$by
    pages <- list()
    failures <- character()
    index$page <- rep(NA_integer_, nrow(index))
    todo <- which(index$status == "ok")
    buffer <- list()

    bar <- if (isTRUE(progress) && length(todo) > 0) {
        utils::txtProgressBar(max = length(todo), style = 3)
    }
    if (!is.null(bar)) {
        on.exit(close(bar), add = TRUE)
    }

    fail <- function(rows, e) {
        failures <<- c(failures, paste0(
            paste(index[[by]][rows], collapse = ", "), ": ",
            conditionMessage(e)))
        index$status[rows] <<- "error"
    }

    flush <- function() {
        if (length(buffer) == 0) {
            return()
        }
        rows <- vapply(buffer, `[[`, integer(1), "row")
        page <- c(
            list(number = length(pages) + 1L, items = buffer),
            .page_tables(buffer, layout))
        buffer <<- list()
        tryCatch({
            page$figure <- render(page)
            # The figure holds what the plots drew; drop them to keep the
            # memory of a large report bounded.
            page$items <- lapply(page$items, function(item) {
                item[setdiff(names(item), c("plot", "ms2_plot"))]
            })
            pages[[page$number]] <<- page
            index$page[rows] <<- page$number
        }, error = function(e) fail(rows, e))
    }

    for (k in seq_along(todo)) {
        i <- todo[k]
        item <- tryCatch(
            .report_item(object, index[i, , drop = FALSE], layout, ms2),
            error = function(e) {
                fail(i, e)
                NULL
            })
        if (!is.null(item)) {
            item$row <- i
            buffer[[length(buffer) + 1]] <- item
            if (length(buffer) == layout$per_page) {
                flush()
            }
        } else if (index$status[i] == "ok") {
            index$status[i] <- "no data"
        }
        if (!is.null(bar)) {
            utils::setTxtProgressBar(bar, k)
        }
    }
    flush()

    if (length(failures) > 0) {
        warning(
            "lp_report: ", length(failures), " page(s) failed to render: ",
            paste(failures, collapse = "; "), call. = FALSE)
    }
    list(index = index, pages = pages)
}

# The plots and data of one feature of the report, or NULL when it has no
# data points, no spectra, and no detected peaks.
.report_item <- function(object, row, layout, ms2) {
    ids <- row$feature_ids[[1]]
    one <- subset_features(object, .data$feature_id %in% !!ids)
    if (.has_no_points(one)) {
        return(NULL)
    }

    ms2_plot <- NULL
    if (!is.null(ms2)) {
        ms2_one <- subset_features(ms2, .data$feature_id %in% !!ids)
        if (length(ms2_one@data@spectra) > 0) {
            ms2_plot <- .render_plot(
                ms2_one, .additional_datasets(ms2_one))@plot
        }
    }

    list(
        key = row[[layout$by]][1],
        heading = .report_heading(row, layout$by),
        plot = .render_plot(one, .additional_datasets(one))@plot,
        ms2_plot = ms2_plot,
        has_ms2 = !is.null(ms2),
        ms2_found = !is.null(ms2_plot),
        annotations = row,
        peaks = .report_peaks(one))
}

.has_no_points <- function(obj) {
    data <- obj@data
    nrow(as_plot_df(data@chromatograms)) == 0 &&
        length(data@spectra) == 0 &&
        nrow(data@detected_peaks) == 0
}
