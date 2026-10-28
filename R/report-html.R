# The figure of an HTML page, as a base64-encoded PNG image. With several
# features, the plots keep the proportions of a single one, with room for
# their titles.
.html_page_figure <- function(page, layout, width, height, dpi) {
    n_items <- length(page$items)
    n_cols <- min(layout$ncol, n_items)
    has_ms2 <- any(vapply(page$items, `[[`, logical(1), "ms2_found"))
    row_height <- if (n_items == 1) height else height / n_cols * 1.15
    if (has_ms2) {
        row_height <- row_height * 5 / 3
    }
    .png_base64(
        .page_figure(page, layout), width,
        row_height * ceiling(n_items / n_cols), dpi)
}

.png_base64 <- function(plot, width, height, dpi) {
    path <- tempfile(fileext = ".png")
    on.exit(unlink(path))
    ggsave(
        path, plot, width = width, height = height, dpi = dpi, units = "in")
    jsonlite::base64_enc(readBin(path, "raw", file.size(path)))
}

.html_table <- function(df, class = NULL, links = NULL) {
    if (is.null(df) || nrow(df) == 0 || ncol(df) == 0) {
        return(NULL)
    }
    rows <- lapply(seq_len(nrow(df)), function(i) {
        cells <- lapply(seq_along(df), function(j) {
            value <- df[[j]][i]
            if (j == 1 && !is.null(links) && !is.na(links[i])) {
                value <- htmltools::tags$a(href = links[i], value)
            } else if (colnames(df)[j] == "status") {
                value <- htmltools::tags$span(
                    class = paste0(
                        "lp-status lp-status-", gsub(" ", "-", value)),
                    value)
            }
            htmltools::tags$td(value)
        })
        htmltools::tags$tr(cells)
    })
    htmltools::div(
        class = "lp-scroll",
        htmltools::tags$table(
            class = class,
            htmltools::tags$thead(htmltools::tags$tr(
                lapply(colnames(df), htmltools::tags$th))),
            htmltools::tags$tbody(rows)))
}

.html_page <- function(page) {
    tags <- htmltools::tags
    heading <- .page_heading(page)
    missing_ms2 <- .missing_ms2(page)
    single <- length(page$items) == 1

    tags$section(
        class = "lp-report-page", id = paste0("page-", page$number),
        tags$h2(heading),
        tags$img(
            src = paste0("data:image/png;base64,", page$figure),
            alt = heading),
        if (length(missing_ms2) > 0) {
            tags$p(class = "lp-note", if (single) {
                "No MS2 spectra match this feature."
            } else {
                paste0(
                    "No MS2 spectra match ",
                    paste(missing_ms2, collapse = ", "), ".")
            })
        },
        if (!is.null(page$annotations) || !is.null(page$peaks)) {
            tags$div(
                class = "lp-tables",
                if (!is.null(page$annotations)) {
                    tags$div(
                        tags$h3("Annotations"),
                        .html_table(
                            page$annotations, if (single) "lp-fields"))
                },
                if (!is.null(page$peaks)) {
                    tags$div(
                        tags$h3("Detected peaks"),
                        if (nrow(page$peaks) > 0) {
                            .html_table(page$peaks)
                        } else {
                            tags$p(class = "lp-note", "No detected peaks.")
                        })
                })
        },
        tags$p(
            class = "lp-back", tags$a(href = "#lp-index", "Back to the index")))
}

# The provenance of an HTML report: tables of the software, the pipeline and
# the inputs, and the full record as JSON for tools to read.
.html_provenance <- function(prov) {
    tags <- htmltools::tags
    tables <- .provenance_tables(prov)
    # "</" would close the script element early.
    json <- gsub("</", "<\\/", .provenance_json(prov), fixed = TRUE)

    tags$section(
        id = "lp-provenance",
        tags$h2("Provenance"),
        tags$p(paste0(
            "The full record, with the version of every loaded package, is ",
            "embedded in this file and written next to it as ",
            prov$output$file, ".provenance.json.")),
        tags$h3("Software"),
        .html_table(tables$software, "lp-fields"),
        tags$h3("Pipeline"),
        .html_table(tables$pipeline),
        if (length(tables$files) > 0) {
            htmltools::tagList(tags$h3("Inputs"), .html_table(tables$files))
        },
        tags$script(
            type = "application/json", id = "lp-provenance-data",
            htmltools::HTML(json)),
        tags$p(
            class = "lp-back", tags$a(href = "#lp-index", "Back to the index")))
}

.write_report_html <- function(rendered, file, title, layout, prov = NULL) {
    tags <- htmltools::tags
    index <- rendered$index
    links <- ifelse(
        is.na(index$page), NA_character_, paste0("#page-", index$page))

    html <- htmltools::tagList(
        tags$head(
            tags$meta(
                name = "viewport",
                content = "width=device-width, initial-scale=1"),
            if (!is.null(prov)) {
                tags$meta(name = "generator", content = prov$generator)
            },
            tags$title(title),
            tags$style(htmltools::HTML(paste0(
                .report_base_css(layout$font_size), layout$style$css)))),
        tags$main(
            tags$h1(title),
            tags$section(
                id = "lp-index",
                tags$h2("Index"),
                tags$p(paste0(
                    nrow(index), " value(s) of '", layout$by, "', ",
                    sum(!is.na(index$page)), " on ",
                    length(rendered$pages), " page(s).")),
                tags$input(
                    id = "lp-search", type = "search",
                    placeholder = "Search the index"),
                .html_table(
                    .report_index_table(index, layout$by, layout$rt_unit),
                    "lp-index-table", links)),
            lapply(rendered$pages, .html_page),
            if (!is.null(prov)) .html_provenance(prov)),
        tags$script(htmltools::HTML(.report_js)))

    htmltools::save_html(html, file)
}

# Filters the index rows on the search box, and sorts them on a click on a
# header (numerically when every value of the column is a number).
.report_js <- "
(function () {
  var table = document.querySelector('.lp-index-table');
  if (!table) return;
  var body = table.tBodies[0];
  document.getElementById('lp-search').addEventListener('input', function () {
    var q = this.value.toLowerCase();
    Array.prototype.forEach.call(body.rows, function (row) {
      row.style.display =
        row.textContent.toLowerCase().indexOf(q) === -1 ? 'none' : '';
    });
  });
  Array.prototype.forEach.call(table.tHead.rows[0].cells, function (th, j) {
    th.addEventListener('click', function () {
      var asc = !th.classList.contains('lp-asc');
      Array.prototype.forEach.call(table.tHead.rows[0].cells, function (h) {
        h.classList.remove('lp-asc', 'lp-desc');
      });
      th.classList.add(asc ? 'lp-asc' : 'lp-desc');
      var rows = Array.prototype.slice.call(body.rows);
      var text = function (r) { return r.cells[j].textContent.trim(); };
      var numeric = rows.every(function (r) {
        return text(r) === '' || !isNaN(Number(text(r)));
      });
      rows.sort(function (a, b) {
        var x = text(a), y = text(b);
        var d = numeric ? (Number(x) || 0) - (Number(y) || 0) :
          x.localeCompare(y);
        return asc ? d : -d;
      });
      rows.forEach(function (r) { body.appendChild(r); });
    });
  });
})();
"
