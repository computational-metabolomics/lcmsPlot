#' UI for the BPC / TIC tab
#'
#' @param id The module namespace id.
#' @return A `shiny::tagList`.
#' @keywords internal
.chromatogram_ui <- function(id) {
    ns <- shiny::NS(id)
    shiny::tags$div(
        class = "lp-grid-content",
        shiny::tags$section(
            class = "lp-plot-card",
            shiny::tags$div(
                class = "flex items-center justify-between mb-4",
                shiny::tags$div(
                    shiny::tags$h2(
                        class = "text-2xl mt-0 mb-2 font-semibold text-slate-800",
                        "Base-peak / total-ion chromatogram"),
                    shiny::tags$p(
                        class = "text-lg text-slate-500",
                        paste("Aggregate each scan and trace the result",
                              "across retention time."))),
                shiny::radioButtons(
                    inputId  = ns("agg"),
                    label    = NULL,
                    choices  = c("BPC (max)" = "max", "TIC (sum)" = "sum"),
                    selected = "max",
                    inline   = TRUE)),
            shiny::plotOutput(ns("plot"), height = "520px")),
        .options_ui(ns("opts"))
    )
}

#' Server for the BPC / TIC tab
#'
#' @param id The module namespace id.
#' @param lcms_obj A reactive returning an `lcmsPlotClass` object, or `NULL`
#'   if no data has been loaded yet.
#' @param selected_samples A reactive returning a `character` vector of
#'   sample IDs to plot, or `NULL` for all samples.
#' @param metadata_cols A reactive returning a `character` vector of metadata
#'   column names to feed the Options panel pickers.
#' @return The module server's return value (currently `NULL`).
#' @keywords internal
.chromatogram_server <- function(id, lcms_obj, selected_samples,
                                 metadata_cols) {
    shiny::moduleServer(id, function(input, output, session) {
        opts <- .options_server("opts", metadata_cols)

        output$plot <- shiny::renderPlot({
            obj <- lcms_obj()
            shiny::req(obj)
            tryCatch({
                layered <- obj +
                    lp_chromatogram(
                        sample_ids = selected_samples(),
                        aggregation_fun = input$agg)
                layered <- .apply_options(layered, opts())
                print(layered + lp_get_plot())
            }, error = function(e) {
                .toast_error(conditionMessage(e),
                             title = "Chromatogram error")
                NULL
            })
        })
    })
}

#' UI for the EIC tab
#'
#' @param id The module namespace id.
#' @return A `shiny::tagList`.
#' @keywords internal
.eic_ui <- function(id) {
    ns <- shiny::NS(id)
    shiny::tags$div(
        class = "lp-grid-content",
        shiny::tags$section(
            class = "lp-plot-card",
            shiny::tags$div(
                class = "mb-4",
                shiny::tags$h2(
                    class = "text-2xl mt-0 mb-2 font-semibold text-slate-800",
                    "Extracted ion chromatogram"),
                shiny::tags$p(
                    class = "text-lg text-slate-500",
                    "Trace intensity for a target m/z window.")),
            shiny::tags$div(
                class = "grid grid-cols-4 gap-3 mb-4",
                shiny::numericInput(ns("mz"),  "m/z",  value = 335,
                                    min = 0, step = 0.001),
                shiny::numericInput(ns("ppm"), "ppm",  value = 20,
                                    min = 0, step = 1),
                shiny::numericInput(ns("rt"),  "RT (sec, optional)",
                                    value = NA, min = 0, step = 1),
                shiny::numericInput(ns("rt_tol"), "RT tolerance (sec)",
                                    value = 30, min = 0, step = 1)),
            shiny::actionButton(ns("go"), "Extract",
                                icon = shiny::icon("crosshairs"),
                                class = "lp-btn-primary mb-4"),
            shiny::plotOutput(ns("plot"), height = "480px")),
        .options_ui(ns("opts"))
    )
}

#' Server for the EIC tab
#'
#' @param id The module namespace id.
#' @param lcms_obj Reactive returning the base `lcmsPlotClass`.
#' @param selected_samples Reactive returning the selected sample IDs.
#' @param metadata_cols Reactive returning available metadata column names.
#' @return `NULL`.
#' @keywords internal
.eic_server <- function(id, lcms_obj, selected_samples, metadata_cols) {
    shiny::moduleServer(id, function(input, output, session) {
        opts <- .options_server("opts", metadata_cols)

        plot_data <- shiny::eventReactive(input$go, {
            obj <- lcms_obj()
            shiny::req(obj)
            shiny::req(input$mz, input$ppm, input$rt_tol)
            tryCatch({
                feats <- .eic_features(
                    mz = input$mz,
                    ppm = input$ppm,
                    rt = if (is.na(input$rt)) NULL else input$rt,
                    rt_tol = input$rt_tol)
                layered <- obj +
                    lp_chromatogram(
                        features    = feats,
                        sample_ids  = selected_samples(),
                        ppm         = input$ppm,
                        rt_tol      = input$rt_tol)
                layered <- .apply_options(layered, opts())
                layered + lp_get_plot()
            }, error = function(e) {
                .toast_error(conditionMessage(e), title = "EIC error")
                NULL
            })
        })

        output$plot <- shiny::renderPlot({
            p <- plot_data()
            shiny::req(p)
            print(p)
        })
    })
}
