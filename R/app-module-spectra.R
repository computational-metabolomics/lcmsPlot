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
                class = "mb-4",
                shiny::tags$h2(
                    class = "text-2xl mt-0 mb-2 font-semibold text-slate-800",
                    "Spectra"),
                shiny::tags$p(
                    class = "text-lg text-slate-500",
                    "Plot a mass spectrum at a given scan index.")),
            shiny::tags$div(
                class = "grid grid-cols-3 gap-3 mb-4",
                shiny::uiOutput(ns("sample")),
                shiny::numericInput(ns("scan_index"), "Scan index",
                                    value = 1, min = 1, step = 1),
                shiny::numericInput(ns("ms_level"),  "MS level",
                                    value = 1, min = 1, step = 1)),
            shiny::actionButton(ns("go"), "Plot",
                                icon = shiny::icon("chart-column"),
                                class = "lp-btn-primary mb-4"),
            shiny::plotOutput(ns("plot"), height = "480px")),
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

        output$sample <- shiny::renderUI({
            obj <- lcms_obj()
            shiny::req(obj)
            ids <- obj@data@metadata$sample_id
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
            p <- plot_data()
            shiny::req(p)
            print(p)
        })
    })
}
