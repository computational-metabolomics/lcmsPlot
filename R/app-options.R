#' Build the per-tab Options panel UI
#'
#' Renders a card with form controls for `lp_facets`, `lp_arrange`,
#' `lp_legend`, and `lp_labels`. Inputs are placed inside the module
#' namespace given by `id`.
#'
#' @param id The module namespace id.
#' @param metadata_cols A `character` vector of metadata column names to
#'   offer as choices for the facet and grouping pickers. Pass `character(0)`
#'   before any data is loaded — the dropdowns will simply be empty.
#' @return A `shiny::tagList` containing the option controls wrapped in a
#'   Tailwind `lp-card`.
#' @keywords internal
.options_ui <- function(id, metadata_cols = character(0)) {
    ns <- shiny::NS(id)
    shiny::tags$aside(
        class = "lp-card-tight",
        shiny::tags$h3(class = "lp-section-title", "Options"),
        shiny::tags$div(
            class = "space-y-3",
            shiny::selectInput(
                ns("facets"),
                "Facet by",
                choices  = metadata_cols,
                selected = NULL,
                multiple = TRUE),
            shiny::numericInput(
                ns("facet_ncol"),
                "Facet columns",
                value = NA, min = 1, step = 1),
            shiny::tags$div(
                class = "flex items-center gap-4",
                shiny::checkboxInput(ns("free_x"), "Free x", FALSE),
                shiny::checkboxInput(ns("free_y"), "Free y", FALSE)),
            shiny::selectInput(
                ns("arrange_by"),
                "Group / colour by",
                choices  = c("(none)" = "", metadata_cols),
                selected = ""),
            shiny::selectInput(
                ns("legend_position"),
                "Legend position",
                choices  = c("default" = "",
                             "top", "right", "bottom", "left", "none"),
                selected = ""),
            shiny::textInput(ns("title"),        "Plot title",  value = ""),
            shiny::textInput(ns("legend_title"), "Legend title", value = "")
        )
    )
}

#' Server-side counterpart for the Options panel
#'
#' Collects the input values into a single reactive list with one entry per
#' option. Empty strings and `NA` are normalised to `NULL` so that
#' [`.apply_options()`] can simply test for non-`NULL` values when deciding
#' which `lp_*()` layer to append.
#'
#' If `metadata_cols` is supplied, the `Facet by` and `Group by` dropdowns
#' are kept in sync with the loaded data's metadata columns.
#'
#' @param id The module namespace id.
#' @param metadata_cols Optional reactive returning a `character` vector of
#'   metadata column names to populate the facet/group dropdowns.
#' @return A reactive returning a named `list` with members `facets`,
#'   `facet_ncol`, `free_x`, `free_y`, `arrange_by`, `legend_position`,
#'   `title`, `legend_title`.
#' @keywords internal
.options_server <- function(id, metadata_cols = NULL) {
    shiny::moduleServer(id, function(input, output, session) {
        if (!is.null(metadata_cols)) {
            shiny::observe({
                cols <- metadata_cols()
                shiny::updateSelectInput(
                    session,
                    inputId  = "facets",
                    choices  = cols,
                    selected = shiny::isolate(input$facets))
                shiny::updateSelectInput(
                    session,
                    inputId  = "arrange_by",
                    choices  = c("(none)" = "", cols),
                    selected = shiny::isolate(input$arrange_by))
            })
        }
        shiny::reactive({
            list(
                facets          = .empty_to_null(input$facets),
                facet_ncol      = .na_to_null(input$facet_ncol),
                free_x          = isTRUE(input$free_x),
                free_y          = isTRUE(input$free_y),
                arrange_by      = .empty_to_null(input$arrange_by),
                legend_position = .empty_to_null(input$legend_position),
                title           = .empty_to_null(input$title),
                legend_title    = .empty_to_null(input$legend_title)
            )
        })
    })
}

#' Apply Options-panel choices as `lp_*()` layers to an lcmsPlotClass
#'
#' Pure function (no Shiny). Given an `lcmsPlotClass` and a list as produced
#' by [`.options_server()`], appends the corresponding `lp_facets`,
#' `lp_arrange`, `lp_legend`, and `lp_labels` layers when their inputs are
#' non-`NULL`. Used by every plot module immediately before `lp_get_plot()`.
#'
#' @param obj An `lcmsPlotClass`.
#' @param opts A named `list` of options (see [`.options_server()`]).
#' @return The (possibly layered) `lcmsPlotClass`.
#' @keywords internal
.apply_options <- function(obj, opts) {
    if (is.null(opts)) return(obj)

    if (!is.null(opts$facets) && length(opts$facets) > 0) {
        obj <- obj + lp_facets(
            facets = opts$facets,
            ncol   = opts$facet_ncol,
            free_x = isTRUE(opts$free_x),
            free_y = isTRUE(opts$free_y))
    }
    if (!is.null(opts$arrange_by)) {
        obj <- obj + lp_arrange(group_by = opts$arrange_by)
    }
    if (!is.null(opts$legend_position)) {
        obj <- obj + lp_legend(position = opts$legend_position)
    }
    if (!is.null(opts$title) || !is.null(opts$legend_title)) {
        obj <- obj + lp_labels(
            title  = opts$title,
            legend = opts$legend_title)
    }
    obj
}

#' Normalise an empty-string input to NULL
#' @keywords internal
.empty_to_null <- function(x) {
    if (is.null(x)) return(NULL)
    if (length(x) == 0) return(NULL)
    if (length(x) == 1 && is.character(x) && !nzchar(x)) return(NULL)
    x
}

#' Normalise an NA numeric input to NULL
#' @keywords internal
.na_to_null <- function(x) {
    if (is.null(x)) return(NULL)
    if (length(x) == 1 && is.na(x)) return(NULL)
    x
}
