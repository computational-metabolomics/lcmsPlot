#' Tabs displayed in the dashboard's sidebar navigation
#'
#' Each entry is the `value` used by the hidden tabsetPanel, the short
#' `label` rendered into the nav links by [`.app_sidebar()`], and the
#' `title` / `description` shown in the topbar while the tab is active.
#'
#' @keywords internal
.APP_TABS <- list(
    list(value       = "chrom",
         label       = "BPC / TIC",
         title       = "Base-peak / total-ion chromatogram",
         description = paste("Aggregate each scan and trace the result",
                             "across retention time.")),
    list(value       = "eic",
         label       = "EIC",
         title       = "Extracted ion chromatogram",
         description = "Trace intensity for a target m/z window."),
    list(value       = "tic",
         label       = "TIC per sample",
         title       = "Total ion current per sample",
         description = paste("Distribution of ion currents per sample.",
                             "Requires XCMSnExp / MsExperiment.")),
    list(value       = "density",
         label       = "Peak density",
         title       = "Peak density",
         description = paste("Density of detected peaks across retention",
                             "time. Requires XCMSnExp / MsExperiment.")),
    list(value       = "spectra",
         label       = "Spectra",
         title       = "Spectra",
         description = "Plot a mass spectrum at a given scan index.")
)

#' Empty-state message shown in plot areas before any data is loaded
#' @keywords internal
.APP_NO_DATA_MSG <- "Upload files and click Load to get started."

#' Empty-state message shown when every sample has been deselected
#' @keywords internal
.APP_NO_SAMPLES_MSG <- "Select at least one sample to plot."

#' Look up an entry of [`.APP_TABS`] by its `value`
#'
#' @param value A `character` tab value, e.g. `"eic"`.
#' @return The matching `list` entry, or `NULL` if there is none.
#' @keywords internal
.app_tab <- function(value) {
    for (t in .APP_TABS) {
        if (identical(t$value, value)) return(t)
    }
    NULL
}

#' Build a titled group for the sidebar
#'
#' @param label A `character` group label, rendered in small caps.
#' @param ... Group content.
#' @param actions Optional tags shown to the right of the label.
#' @param inset If `TRUE`, the content is indented to line up with the
#'   label and the nav link text.
#' @return A `shiny::tags$section`.
#' @keywords internal
.sidebar_group <- function(label, ..., actions = NULL, inset = TRUE) {
    header <- if (is.null(actions)) {
        shiny::tags$h2(class = "lp-group-label", label)
    } else {
        shiny::tags$div(
            class = "lp-group-header",
            shiny::tags$h2(class = "lp-group-label", label),
            shiny::tags$div(class = "flex items-center gap-3", actions))
    }
    body <- shiny::tags$div(
        class = if (inset) "lp-sidebar-body" else "flex flex-col gap-1",
        ...)
    shiny::tags$section(class = "lp-sidebar-group", header, body)
}

#' Build the sticky topbar of the dashboard
#'
#' Shows the active view's title and description (rendered server-side
#' into `page_header`) and a chip with the number of loaded samples.
#'
#' @return A `shiny::tags$header`.
#' @keywords internal
.app_topbar <- function() {
    shiny::tags$header(
        class = "lp-topbar",
        shiny::uiOutput("page_header"),
        shiny::uiOutput("samples_chip", inline = TRUE)
    )
}

#' Build the sidebar (wordmark, view navigation, file upload, samples)
#'
#' @return A `shiny::tags$aside`.
#' @keywords internal
.app_sidebar <- function() {
    shiny::tags$aside(
        class = "lp-sidebar",
        shiny::tags$p(
            class = "lp-wordmark",
            shiny::HTML('lcms<span class="text-brand-500">Plot</span>')),
        shiny::tags$nav(
            `aria-label` = "Views",
            .sidebar_group(
                "Views",
                inset = FALSE,
                lapply(seq_along(.APP_TABS), function(i) {
                    t <- .APP_TABS[[i]]
                    link <- shiny::actionLink(
                        inputId = paste0("nav_", t$value),
                        label   = t$label,
                        class   = if (i == 1L) {
                            "lp-nav-item lp-nav-active"
                        } else {
                            "lp-nav-item"
                        })
                    if (i == 1L) {
                        link <- shiny::tagAppendAttributes(
                            link, `aria-current` = "page")
                    }
                    link
                }))),
        .sidebar_group(
            "Data",
            shiny::fileInput(
                inputId = "files",
                label   = "Files",
                multiple = TRUE,
                accept  = c(".mzML", ".mzml", ".CDF", ".cdf",
                            ".raw", ".RAW",
                            ".cdResult",
                            ".rds", ".Rds",
                            ".RData", ".Rdata", ".rdata", ".rda")),
            shiny::tags$p(
                class = "lp-help",
                paste("One or more mzML, CDF or raw files, or a single",
                      "cdResult, rds or RData file.")),
            shiny::actionButton(
                inputId = "load",
                label   = "Load",
                class   = "lp-btn-primary w-full")
        ),
        .sidebar_group(
            "Samples",
            actions = shiny::tagList(
                shiny::actionLink("samples_all",  "All",
                                  class = "lp-link-btn"),
                shiny::actionLink("samples_none", "None",
                                  class = "lp-link-btn")),
            shiny::tags$div(
                class = "lp-sample-list",
                shiny::uiOutput("sample_picker"))
        )
    )
}

#' Build the content area: hidden tabset of plot views
#'
#' Each [shiny::tabPanel()] hosts a plot card on the left and the Options
#' card on the right. Tab switching is driven from the sidebar navigation
#' via [shiny::updateTabsetPanel()] — the tabsetPanel itself is
#' `type = "hidden"` so no Bootstrap tab strip is shown.
#'
#' @return A `shiny::tags$main`.
#' @keywords internal
.app_content <- function() {
    shiny::tags$main(
        class = "lp-content",
        shiny::tabsetPanel(
            id = "tabs",
            type = "hidden",
            shiny::tabPanel(
                "chrom",
                .chromatogram_ui("chrom")),
            shiny::tabPanel(
                "eic",
                .eic_ui("eic")),
            shiny::tabPanel(
                "tic",
                .tic_ui("tic")),
            shiny::tabPanel(
                "density",
                .peak_density_ui("density")),
            shiny::tabPanel(
                "spectra",
                .spectra_ui("spectra"))
        )
    )
}
