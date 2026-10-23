#' Keep the data matching a condition on its annotations
#'
#' `subset_features()` narrows the data already extracted for an
#' `lcmsPlotClass` object, without reading the raw data again. Each condition
#' is applied to every dataset that carries at least one of the columns it
#' refers to: the chromatograms (through their `chromData()` columns), the
#' spectra (through their spectra variables), the detected peaks, and the
#' summaries. Datasets that carry none of them are kept whole.
#'
#' This is useful to draw one plot per compound or feature out of a single
#' extraction.
#'
#' @param obj An instance of class `lcmsPlotClass`.
#' @param ... Conditions on the dataset columns, as in [dplyr::filter()].
#' @return The `lcmsPlotClass` object with the matching data only.
#' @export
#' @examples
#' raw_files <- dir(
#'    system.file("cdf", package = "faahKO"),
#'    full.names = TRUE,
#'    recursive = TRUE)[1:2]
#'
#' features <- data.frame(
#'     feature_id = c("A", "B"),
#'     mzmin = c(334.9, 344.0),
#'     mzmax = c(335.1, 344.2),
#'     rtmin = c(2700, 3400),
#'     rtmax = c(2900, 3500))
#'
#' p <- lcmsPlot(raw_files) + lp_chromatogram(features = features)
#' p_a <- subset_features(p, feature_id == "A")
subset_features <- function(obj, ...) {
    conditions <- rlang::enquos(...)
    referenced <- unique(unlist(lapply(
        conditions, function(q) all.vars(rlang::quo_get_expr(q)))))

    matching_rows <- function(df) {
        if (nrow(df) == 0 || !any(referenced %in% colnames(df))) {
            return(NULL)
        }
        df |>
            as_tibble() |>
            mutate(.row = row_number()) |>
            filter(!!!conditions) |>
            pull(".row")
    }

    data <- obj@data

    if (length(data@chromatograms) > 0) {
        rows <- matching_rows(as.data.frame(
            Chromatograms::chromData(data@chromatograms)))
        if (!is.null(rows)) {
            data@chromatograms <- data@chromatograms[rows]
        }
    }

    if (length(data@spectra) > 0) {
        variables <- intersect(
            referenced, Spectra::spectraVariables(data@spectra))
        if (length(variables) > 0) {
            rows <- matching_rows(as.data.frame(
                Spectra::spectraData(data@spectra, variables)))
            data@spectra <- data@spectra[rows]
        }
    }

    rows <- matching_rows(data@detected_peaks)
    if (!is.null(rows)) {
        data@detected_peaks <- data@detected_peaks[rows, , drop = FALSE]
    }

    data@summaries <- lapply(data@summaries, function(df) {
        rows <- matching_rows(df)
        if (is.null(rows)) df else df[rows, , drop = FALSE]
    })

    validObject(data)
    obj@data <- data
    obj
}
