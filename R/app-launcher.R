#' `htmlDependency` for the bundled Tailwind stylesheet
#'
#' Standard Shiny pattern for serving a CSS file from `inst/www/`. Uses
#' [htmltools::htmlDependency()] so the stylesheet picks up a version stamp
#' (bust the browser cache on package upgrades) and resolves the path via
#' [system.file()] internally — works under [devtools::load_all()] too.
#'
#' @keywords internal
.lcmsPlot_dep <- function() {
    htmltools::htmlDependency(
        name       = "lcmsPlot",
        version    = as.character(utils::packageVersion("lcmsPlot")),
        package    = "lcmsPlot",
        src        = "www",
        stylesheet = "lcmsPlot.css")
}

#' Internal UI builder for the Shiny app
#'
#' Lays out the dashboard: sticky topbar, left nav rail, main content with
#' a hidden tabset. The Tailwind stylesheet is attached via
#' [`.lcmsPlot_dep()`].
#'
#' @return A `shiny::tagList` wrapping the dashboard.
#' @keywords internal
.build_ui <- function() {
    shiny::tagList(
        .lcmsPlot_dep(),
        shiny::tags$head(
            shiny::tags$meta(charset = "utf-8"),
            shiny::tags$meta(
                name    = "viewport",
                content = "width=device-width, initial-scale=1"),
            shiny::tags$title("lcmsPlot"),
            if (requireNamespace("shinytoastr", quietly = TRUE)) {
                shinytoastr::useToastr()
            }
        ),
        shiny::tags$div(
            class = "min-h-screen bg-slate-50 text-slate-800",
            .app_topbar(),
            shiny::tags$div(
                class = "lp-grid-main",
                .app_navrail(),
                .app_content()
            )
        )
    )
}

#' Internal server builder for the Shiny app
#'
#' @return A function suitable for `shiny::shinyApp(server = ...)`.
#' @keywords internal
.build_server <- function() {
    function(input, output, session) {
        # nocov start
        session_dir <- file.path(tempdir(), paste0("lcmsPlot-", session$token))

        lcms_obj <- shiny::eventReactive(input$load, {
            paths <- .uploaded_files_to_paths(input$files, session_dir)
            if (length(paths) == 0) {
                .toast_error("Upload at least one file before loading.",
                             title = "No files")
                return(NULL)
            }
            tryCatch({
                dataset <- .load_dataset(paths)
                lcmsPlot(dataset)
            }, error = function(e) {
                .toast_error(conditionMessage(e), title = "Load error")
                NULL
            })
        })

        metadata_cols <- shiny::reactive({
            obj <- lcms_obj()
            if (is.null(obj)) character(0)
            else colnames(obj@data@metadata)
        })

        output$sample_picker <- shiny::renderUI({
            obj <- lcms_obj()
            if (is.null(obj)) {
                return(shiny::tags$p(
                    class = "text-xs text-slate-400",
                    "Load files to see samples."))
            }
            ids <- obj@data@metadata$sample_id
            shiny::checkboxGroupInput(
                inputId  = "samples",
                label    = NULL,
                choices  = ids,
                selected = ids)
        })

        output$samples_chip <- shiny::renderUI({
            obj <- lcms_obj()
            n <- if (is.null(obj)) 0L else length(obj@data@metadata$sample_id)
            shiny::tags$span(
                class = "lp-chip",
                shiny::icon("vial"),
                sprintf("%d sample%s loaded", n, if (n == 1) "" else "s"))
        })

        selected_samples <- shiny::reactive({
            sel <- input$samples
            if (is.null(sel) || length(sel) == 0) NULL else sel
        })

        # Nav-rail observers: switch the hidden tabset when a nav link is
        # clicked, and toggle the `lp-nav-active` class on the link itself.
        active_tab <- shiny::reactiveVal("chrom")
        lapply(.APP_TABS, function(t) {
            shiny::observeEvent(input[[paste0("nav_", t$value)]], {
                shiny::updateTabsetPanel(session, "tabs",
                                         selected = t$value)
                active_tab(t$value)
            })
        })
        shiny::observe({
            current <- active_tab()
            for (t in .APP_TABS) {
                shinyjs <- session$sendCustomMessage
                # Manipulate the active class via a tiny inline JS message.
                session$sendCustomMessage(
                    type = "lp_set_nav_active",
                    message = list(id = paste0("nav_", t$value),
                                   active = identical(t$value, current)))
            }
        })

        # Modules: each gets the base obj, the selected samples, and the
        # metadata columns to populate its Options panel.
        .chromatogram_server("chrom",   lcms_obj, selected_samples,
                             metadata_cols)
        .eic_server         ("eic",     lcms_obj, selected_samples,
                             metadata_cols)
        .tic_server         ("tic",     lcms_obj, selected_samples,
                             metadata_cols)
        .peak_density_server("density", lcms_obj, metadata_cols)
        .spectra_server     ("spectra", lcms_obj, metadata_cols)

        session$onSessionEnded(function() {
            if (dir.exists(session_dir)) unlink(session_dir, recursive = TRUE)
        })
        # nocov end
    }
}

#' Tiny JS handler that toggles the `lp-nav-active` class on a nav link
#'
#' Injected once into the page head so the R `.build_server()` can switch
#' the highlighted nav button via `session$sendCustomMessage()`.
#'
#' @keywords internal
.nav_active_js <- "
Shiny.addCustomMessageHandler('lp_set_nav_active', function(m) {
    var el = document.getElementById(m.id);
    if (!el) return;
    if (m.active) { el.classList.add('lp-nav-active'); }
    else          { el.classList.remove('lp-nav-active'); }
});
"

#' Launch the lcmsPlot interactive Shiny app
#'
#' Returns a [shiny::shinyApp] object that lets users upload one or more raw
#' files (mzML / CDF) and explore them through the existing `lp_*()` plotting
#' API: base-peak / total-ion chromatograms, extracted ion chromatograms,
#' peak density, and spectra by scan index. Each plot tab carries a side
#' Options panel for `lp_facets`, `lp_arrange`, `lp_legend`, and `lp_labels`.
#'
#' Per Bioconductor's Shiny guidelines the function does *not* call
#' [shiny::runApp()] itself — callers run the returned app with
#' `shiny::runApp(lcmsPlotApp())` or by printing it at the R prompt.
#'
#' Requires the optional packages `shiny` and `shinytoastr`; both are
#' declared in `Suggests`. An informative error is raised if either is
#' missing.
#'
#' @param max_upload_size A `numeric` value in bytes setting the maximum
#'   per-file upload size. Defaults to `2 * 1024^3` (2 GB) since raw mzML
#'   files routinely exceed Shiny's 5 MB default. Set with care if running
#'   the app on shared infrastructure.
#' @return A [shiny::shinyApp] object.
#' @export
#' @examples
#' if (interactive()) {
#'     shiny::runApp(lcmsPlotApp())
#' }
lcmsPlotApp <- function(max_upload_size = 2 * 1024^3) {
    if (!requireNamespace("shiny", quietly = TRUE)) {
        stop("`shiny` is required to run lcmsPlotApp(). ",
             "Install it with: install.packages('shiny')",
             call. = FALSE)
    }
    if (!requireNamespace("shinytoastr", quietly = TRUE)) {
        stop("`shinytoastr` is required to run lcmsPlotApp(). ",
             "Install it with: install.packages('shinytoastr')",
             call. = FALSE)
    }
    options(shiny.maxRequestSize = max_upload_size)

    ui <- shiny::bootstrapPage(
        shiny::tags$head(shiny::tags$script(shiny::HTML(.nav_active_js))),
        .build_ui())

    shiny::shinyApp(ui = ui, server = .build_server())
}
