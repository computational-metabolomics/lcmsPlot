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
                class = "flex items-center justify-between mb-4",
                shiny::tags$div(
                    shiny::tags$h2(
                        class = "text-2xl mt-0 mb-2 font-semibold text-slate-800",
                        "Total ion current per sample"),
                    shiny::tags$p(
                        class = "text-lg text-slate-500",
                        paste("Distribution of ion currents per sample.",
                              "Requires XCMSnExp / MsExperiment."))),
                shiny::radioButtons(
                    inputId  = ns("type"),
                    label    = NULL,
                    choices  = c("Boxplot" = "boxplot",
                                 "Violin"  = "violin",
                                 "Jitter"  = "jitter"),
                    selected = "boxplot",
                    inline   = TRUE)),
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
            shiny::req(obj)
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
        })
    })
}
