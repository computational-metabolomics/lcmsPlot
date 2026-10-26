#' UI for the spectra tab
#'
#' @param id The module namespace id.
#' @return A `shiny::tagList`.
#' @keywords internal
.spectra_ui <- function(id) {
    ns <- shiny::NS(id)
    shiny::tags$div(
        class = "lp-grid-content",
        shiny::tags$section(
            class = "lp-plot-card",
            shiny::tags$div(
                class = "lp-toolbar",
                shiny::uiOutput(ns("sample")),
                shiny::numericInput(ns("scan_index"), "Scan index",
                                    value = 1, min = 1, step = 1),
                shiny::numericInput(ns("ms_level"),  "MS level",
                                    value = 1, min = 1, step = 1),
                shiny::actionButton(ns("go"), "Plot",
                                    class = "lp-btn-primary")),
            shiny::plotOutput(ns("plot"), height = "520px")),
        .options_ui(ns("opts"))
    )
}

#' Server for the spectra tab
#'
#' @param id The module namespace id.
#' @param lcms_obj Reactive returning the base `lcmsPlotClass`.
#' @param metadata_cols Reactive returning available metadata column names.
#' @return `NULL`.
#' @keywords internal
.spectra_server <- function(id, lcms_obj, metadata_cols) {
    shiny::moduleServer(id, function(input, output, session) {
        opts <- .options_server("opts", metadata_cols)

        # Rendered before data is loaded too (with no choices) so the
        # toolbar layout stays stable.
        output$sample <- shiny::renderUI({
            obj <- lcms_obj()
            ids <- if (is.null(obj)) character(0)
                   else sampleMetadata(obj)$sample_id
            shiny::selectInput(
                inputId  = session$ns("sample_id"),
                label    = "Sample",
                choices  = ids,
                selected = ids[1])
        })

        plot_data <- shiny::eventReactive(input$go, {
            obj <- lcms_obj()
            shiny::req(obj, input$sample_id, input$scan_index)
            tryCatch({
                layered <- obj +
                    lp_spectra(
                        sample_ids = input$sample_id,
                        scan_index = as.integer(input$scan_index),
                        ms_level   = as.integer(input$ms_level))
                layered <- .apply_options(layered, opts())
                layered + lp_get_plot()
            }, error = function(e) {
                .toast_error(conditionMessage(e), title = "Spectra error")
                NULL
            })
        })

        output$plot <- shiny::renderPlot({
            shiny::validate(shiny::need(lcms_obj(), .APP_NO_DATA_MSG))
            shiny::validate(shiny::need(
                input$go > 0, "Pick a sample and scan, then click Plot."))
            p <- plot_data()
            shiny::req(p)
            print(p)
        }, res = 96)
    })
}
