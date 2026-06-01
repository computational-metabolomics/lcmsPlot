#' Tabs displayed in the dashboard's left nav rail
#'
#' Each entry is the `value` used by the hidden tabsetPanel; the label and
#' optional FontAwesome icon (via [shiny::icon()]) are rendered into the
#' nav buttons by [`.app_navrail()`].
#'
#' @keywords internal
.APP_TABS <- list(
    list(value = "chrom",   label = "BPC / TIC",         icon = "wave-square"),
    list(value = "eic",     label = "EIC",               icon = "crosshairs"),
    list(value = "tic",     label = "TIC per sample",    icon = "chart-simple"),
    list(value = "density", label = "Peak density",      icon = "chart-area"),
    list(value = "spectra", label = "Spectra",           icon = "chart-column")
)

#' Build the sticky topbar of the dashboard
#'
#' @return A `shiny::tags$header`.
#' @keywords internal
.app_topbar <- function() {
    shiny::tags$header(
        class = paste(
            "sticky top-0 z-10 bg-white border-b border-slate-200",
            "px-5 py-3 flex items-center justify-between"),
        shiny::tags$div(
            class = "flex items-center gap-3",
            shiny::tags$div(
                shiny::tags$h1(
                    class = "text-2xl m-0 font-semibold text-slate-800",
                    "lcmsPlot"),
                shiny::tags$p(
                    class = "text-lg m-0 text-slate-500",
                    "Interactive LC-MS data explorer"))
        ),
        shiny::tags$div(
            class = "flex items-center gap-3",
            shiny::uiOutput("samples_chip", inline = TRUE)
        )
    )
}

#' Build the left nav rail (file upload, samples, tab nav)
#'
#' @return A `shiny::tags$aside`.
#' @keywords internal
.app_navrail <- function() {
    shiny::tags$aside(
        class = "flex flex-col gap-6",
        shiny::tags$section(
            class = "lp-card",
            shiny::tags$h2(class = "lp-section-title", "Data"),
            shiny::fileInput(
                inputId = "files",
                label   = "Upload data",
                multiple = TRUE,
                accept  = c(".mzML", ".mzml", ".CDF", ".cdf",
                            ".raw", ".RAW",
                            ".cdResult",
                            ".rds", ".Rds",
                            ".RData", ".Rdata", ".rdata", ".rda")),
            shiny::tags$p(
                class = "text-sm text-slate-500 mt-1 mb-3",
                paste("mzML, CDF, raw (multi-file),",
                      "or a single cdResult / rds / RData.")),
            shiny::actionButton(
                inputId = "load",
                label   = "Load",
                icon    = shiny::icon("upload"),
                class   = "lp-btn-primary w-full")
        ),
        shiny::tags$section(
            class = "lp-card",
            shiny::tags$h2(class = "lp-section-title", "Samples"),
            shiny::uiOutput("sample_picker")
        ),
        shiny::tags$nav(
            class = "lp-card flex flex-col gap-1",
            shiny::tags$h2(class = "lp-section-title", "Views"),
            lapply(seq_along(.APP_TABS), function(i) {
                t <- .APP_TABS[[i]]
                cls <- if (i == 1L) "lp-nav-item lp-nav-active"
                       else        "lp-nav-item"
                shiny::actionLink(
                    inputId = paste0("nav_", t$value),
                    label   = shiny::tagList(
                        shiny::icon(t$icon),
                        shiny::tags$span(t$label)),
                    class   = cls)
            })
        )
    )
}

#' Build the content area: hidden tabset + per-tab options panel
#'
#' Each [shiny::tabPanel()] hosts a plot card on the left and the Options
#' card on the right. Tab switching is driven from the left nav rail via
#' [shiny::updateTabsetPanel()] — the tabsetPanel itself is `type = "hidden"`
#' so no Bootstrap tab strip is shown.
#'
#' @return A `shiny::tags$main`.
#' @keywords internal
.app_content <- function() {
    shiny::tags$main(
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
