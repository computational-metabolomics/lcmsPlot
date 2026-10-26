#' UI for the TIC summary tab
#'
#' @param id The module namespace id.
#' @return A `shiny::tagList`.
#' @keywords internal
.tic_ui <- function(id) {
    ns <- shiny::NS(id)
    shiny::tags$div(
        class = "lp-grid-content",
        shiny::tags$section(
            class = "lp-plot-card",
            shiny::tags$div(
                class = "lp-toolbar",
                shiny::tags$div(
                    class = "lp-segmented",
                    shiny::radioButtons(
                        inputId  = ns("type"),
                        label    = "Plot type",
                        choices  = c("Boxplot" = "boxplot",
                                     "Violin"  = "violin",
                                     "Jitter"  = "jitter"),
                        selected = "boxplot",
                        inline   = TRUE))),
            shiny::plotOutput(ns("plot"), height = "520px")),
        .options_ui(ns("opts"))
    )
}

#' Server for the TIC summary tab
#'
#' @param id The module namespace id.
#' @param lcms_obj Reactive returning the base `lcmsPlotClass`.
#' @param selected_samples Reactive returning the selected sample IDs.
#' @param metadata_cols Reactive returning available metadata column names.
#' @return `NULL`.
#' @keywords internal
.tic_server <- function(id, lcms_obj, selected_samples, metadata_cols) {
    shiny::moduleServer(id, function(input, output, session) {
        opts <- .options_server("opts", metadata_cols)

        output$plot <- shiny::renderPlot({
            obj <- lcms_obj()
            shiny::validate(shiny::need(obj, .APP_NO_DATA_MSG))
            shiny::validate(shiny::need(length(selected_samples()) > 0,
                                        .APP_NO_SAMPLES_MSG))
            tryCatch({
                layered <- obj +
                    lp_total_ion_current(
                        sample_ids = selected_samples(),
                        type       = input$type)
                layered <- .apply_options(layered, opts())
                print(layered + lp_get_plot())
            }, error = function(e) {
                .toast_error(conditionMessage(e), title = "TIC error")
                NULL
            })
        }, res = 96)
    })
}
