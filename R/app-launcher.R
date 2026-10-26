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
#' Lays out the dashboard: a sticky sidebar (navigation, data, samples)
#' next to the main column (sticky topbar, hidden tabset). The Tailwind
#' stylesheet is attached separately by [lcmsPlotApp()]; Shiny's busy indicators
#' add a spinner to recalculating outputs and a progress pulse at the top
#' of the page whenever the server is busy.
#'
#' @return A `shiny::tagList` wrapping the dashboard.
#' @keywords internal
.build_ui <- function() {
    shiny::tagList(
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
        # Top-level so the spinner options scope to <body>.
        shiny::useBusyIndicators(spinners = TRUE, pulse = TRUE),
        shiny::busyIndicatorOptions(
            spinner_type     = "ring",
            spinner_color    = "#465fff",
            spinner_size     = "36px",
            pulse_background = "linear-gradient(90deg, #9cb9ff, #465fff)",
            pulse_height     = "3px"),
        shiny::tags$div(
            class = "lp-shell",
            .app_sidebar(),
            shiny::tags$div(
                class = "lp-main",
                .app_topbar(),
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

        # A reactiveVal rather than an eventReactive so that, before the
        # first load, readers get `NULL` (and can show an empty state)
        # instead of a silent `req()` failure.
        lcms_obj <- shiny::reactiveVal(NULL)
        shiny::observeEvent(input$load, {
            paths <- .uploaded_files_to_paths(input$files, session_dir)
            if (length(paths) == 0) {
                .toast_error("Upload at least one file before loading.",
                             title = "No files")
                lcms_obj(NULL)
                return()
            }
            lcms_obj(tryCatch({
                dataset <- .load_dataset(paths)
                obj <- lcmsPlot(dataset)
                n <- length(sampleMetadata(obj)$sample_id)
                .toast_success(
                    sprintf("Loaded %d sample%s.", n, if (n == 1) "" else "s"),
                    title = "Data loaded")
                obj
            }, error = function(e) {
                .toast_error(conditionMessage(e), title = "Load error")
                NULL
            }))
        })

        sample_ids <- shiny::reactive({
            obj <- lcms_obj()
            if (is.null(obj)) character(0) else sampleMetadata(obj)$sample_id
        })

        metadata_cols <- shiny::reactive({
            obj <- lcms_obj()
            if (is.null(obj)) character(0)
            else colnames(sampleMetadata(obj))
        })

        output$sample_picker <- shiny::renderUI({
            ids <- sample_ids()
            if (length(ids) == 0) {
                return(shiny::tags$p(
                    class = "lp-help",
                    "Load files to see samples."))
            }
            shiny::checkboxGroupInput(
                inputId  = "samples",
                label    = NULL,
                choices  = ids,
                selected = ids)
        })

        output$samples_chip <- shiny::renderUI({
            n <- length(sample_ids())
            if (n == 0) {
                return(shiny::tags$span(class = "lp-chip-muted",
                                        "No data loaded"))
            }
            shiny::tags$span(
                class = "lp-chip",
                sprintf("%d sample%s loaded", n, if (n == 1) "" else "s"))
        })

        shiny::observeEvent(input$samples_all, {
            shiny::updateCheckboxGroupInput(session, "samples",
                                            selected = sample_ids())
        })
        shiny::observeEvent(input$samples_none, {
            shiny::updateCheckboxGroupInput(session, "samples",
                                            selected = character(0))
        })

        # `NULL` until data is loaded; afterwards the checked sample IDs,
        # which may be empty (the modules then show an empty state). IDs
        # left over from a previously loaded dataset are dropped.
        selected_samples <- shiny::reactive({
            ids <- sample_ids()
            if (length(ids) == 0) return(NULL)
            intersect(input$samples, ids)
        })

        # Sidebar nav observers: switch the hidden tabset when a nav link
        # is clicked, and toggle the `lp-nav-active` class on the link.
        active_tab <- shiny::reactiveVal("chrom")
        lapply(.APP_TABS, function(t) {
            shiny::observeEvent(input[[paste0("nav_", t$value)]], {
                shiny::updateTabsetPanel(session, "tabs",
                                         selected = t$value)
                active_tab(t$value)
            })
        })

        output$page_header <- shiny::renderUI({
            t <- .app_tab(active_tab())
            shiny::tags$div(
                shiny::tags$h1(class = "lp-page-title", t$title),
                shiny::tags$p(class = "lp-page-subtitle", t$description))
        })

        shiny::observe({
            current <- active_tab()
            for (t in .APP_TABS) {
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

#' Tiny JS handlers for nav state and button busy state
#'
#' Injected once into the page head. `lp_set_nav_active` lets the R
#' `.build_server()` switch the highlighted nav link (class plus
#' `aria-current`) via `session$sendCustomMessage()`. The click handler
#' marks a primary action button as busy (spinner, disabled) until the
#' server reports it is idle again.
#'
#' @keywords internal
.nav_active_js <- "
Shiny.addCustomMessageHandler('lp_set_nav_active', function(m) {
    var el = document.getElementById(m.id);
    if (!el) return;
    if (m.active) {
        el.classList.add('lp-nav-active');
        el.setAttribute('aria-current', 'page');
    } else {
        el.classList.remove('lp-nav-active');
        el.removeAttribute('aria-current');
    }
});
document.addEventListener('click', function(e) {
    var btn = e.target.closest('.btn.lp-btn-primary.action-button');
    if (!btn) return;
    // Defer so Shiny registers the click before the button is disabled.
    setTimeout(function() {
        btn.classList.add('is-busy');
        btn.setAttribute('aria-busy', 'true');
        btn.disabled = true;
    }, 0);
});
$(document).on('shiny:idle', function() {
    document.querySelectorAll('.btn.is-busy').forEach(function(btn) {
        btn.classList.remove('is-busy');
        btn.removeAttribute('aria-busy');
        btn.disabled = false;
    });
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

    # The stylesheet is a sibling *after* the page so that it loads after
    # Bootstrap and selectize and wins over their equal-specificity rules.
    ui <- shiny::tagList(
        shiny::bootstrapPage(
            shiny::tags$head(shiny::tags$script(shiny::HTML(.nav_active_js))),
            .build_ui()),
        .lcmsPlot_dep())

    shiny::shinyApp(ui = ui, server = .build_server())
}
