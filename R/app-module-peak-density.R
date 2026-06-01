#' UI for the peak density tab
#'
#' @param id The module namespace id.
#' @return A `shiny::tagList`.
#' @keywords internal
.peak_density_ui <- function(id) {
    ns <- shiny::NS(id)
    shiny::tags$div(
        class = "lp-grid-content",
        shiny::tags$section(
            class = "lp-plot-card",
            shiny::tags$div(
                class = "mb-4",
                shiny::tags$h2(
                    class = "text-2xl mt-0 mb-2 font-semibold text-slate-800",
                    "Peak density"),
                shiny::tags$p(
                    class = "text-lg text-slate-500",
                    paste("Requires data with detected peaks",
                          "(XCMSnExp / MsExperiment)."))),
            shiny::tags$div(
                class = "grid grid-cols-3 gap-3 mb-4",
                shiny::numericInput(ns("mz"),  "m/z (optional)",
                                    value = NA,  min = 0, step = 0.001),
                shiny::numericInput(ns("ppm"), "ppm",
                                    value = 20,  min = 0, step = 1),
                shiny::numericInput(ns("bw"),  "Bandwidth (sec)",
                                    value = 30,  min = 1, step = 1)),
            shiny::actionButton(ns("go"), "Plot",
                                icon = shiny::icon("chart-area"),
                                class = "lp-btn-primary mb-4"),
            shiny::plotOutput(ns("plot"), height = "480px")),
        .options_ui(ns("opts"))
    )
}

#' Server for the peak density tab
#'
#' @param id The module namespace id.
#' @param lcms_obj Reactive returning the base `lcmsPlotClass`.
#' @param metadata_cols Reactive returning available metadata column names.
#' @return `NULL`.
#' @keywords internal
.peak_density_server <- function(id, lcms_obj, metadata_cols) {
    shiny::moduleServer(id, function(input, output, session) {
        opts <- .options_server("opts", metadata_cols)

        plot_data <- shiny::eventReactive(input$go, {
            obj <- lcms_obj()
            shiny::req(obj)
            shiny::req(input$bw)
            tryCatch({
                feats <- if (is.na(input$mz)) {
                    NULL
                } else {
                    .eic_features(
                        mz = input$mz,
                        ppm = input$ppm,
                        rt = NULL)
                }
                layered <- obj +
                    lp_peak_density(features = feats, bw = input$bw)
                layered <- .apply_options(layered, opts())
                layered + lp_get_plot()
            }, error = function(e) {
                .toast_error(conditionMessage(e),
                             title = "Peak density error")
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
