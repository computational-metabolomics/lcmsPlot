#' The index of a report
#'
#' One row per value of `by`, collapsed from the chromatograms' annotations:
#' a column keeps its value when the features of the row share one, and lists
#' the distinct values otherwise.
#'
#' @param object An instance of class `lcmsPlotClass`.
#' @param by The annotation column that defines a page.
#' @param features The values of `by` to keep, or `NULL` for all.
#' @return A `tibble` with the `by` column (as `character`), the annotation
#' columns, `feature_ids` (a list-column), `n_samples`, and `status`.
#' @keywords internal
.report_index <- function(object, by, features = NULL) {
    chrom_data <- as.data.frame(
        Chromatograms::chromData(object@data@chromatograms))
    if (!"feature_id" %in% colnames(chrom_data)) {
        stop(
            "lp_report: the data has no features; plot the chromatograms ",
            "of features, e.g. with lp_chromatogram(features = ).",
            call. = FALSE)
    }

    annotations <- setdiff(
        .chrom_annotation_columns(chrom_data), c("sample_index", "rt_type"))
    if (!by %in% annotations) {
        stop(
            "lp_report: 'by' must be one of the annotation columns: ",
            paste0("'", annotations, "'", collapse = ", "), ".",
            call. = FALSE)
    }

    chrom_data <- chrom_data[!is.na(chrom_data[[by]]), , drop = FALSE]
    keys <- as.character(chrom_data[[by]])
    groups <- factor(keys, levels = unique(keys))
    index <- tibble(!!by := levels(groups))
    index$feature_ids <- unname(lapply(
        split(as.character(chrom_data$feature_id), groups), unique))
    for (column in setdiff(annotations, c(by, "feature_id"))) {
        index[[column]] <- .collapse_by_group(chrom_data[[column]], groups)
    }

    index <- .add_peak_summaries(index, object@data@detected_peaks)
    # Without peaks to locate them, the features sit at the middle of their
    # extraction windows.
    windows <- list(mz = c("mzMin", "mzMax"), rt = c("rtMin", "rtMax"))
    for (column in names(windows)) {
        bounds <- windows[[column]]
        if (!column %in% colnames(index) &&
            all(bounds %in% colnames(chrom_data))) {
            middle <- rowMeans(chrom_data[, bounds, drop = FALSE])
            index[[column]] <- unname(vapply(
                split(middle, groups), stats::median, numeric(1)))
        }
    }
    index$status <- rep("ok", nrow(index))

    if (!is.null(features)) {
        features <- unique(as.character(features))
        missing <- setdiff(features, index[[by]])
        if (length(missing) > 0) {
            index <- bind_rows(index, tibble(
                !!by := missing,
                feature_ids = rep(list(character()), length(missing)),
                n_samples = 0L,
                status = "no data"))
        }
        index <- index[match(features, index[[by]]), , drop = FALSE]
    }

    first <- intersect(
        c(by, "name", "formula", "adduct", "adduct_ion", "mz", "rt"),
        colnames(index))
    index[, c(first, setdiff(colnames(index), first))]
}

# One value per group: the shared value, or the distinct values pasted
# together when a group holds several.
.collapse_by_group <- function(x, groups) {
    values <- lapply(split(x, groups), function(v) unique(v[!is.na(v)]))
    if (all(lengths(values) <= 1)) {
        return(unname(do.call(c, lapply(values, function(v) {
            if (length(v) == 0) x[NA_integer_] else v
        }))))
    }
    unname(vapply(values, function(v) {
        paste(sort(as.character(v)), collapse = "; ")
    }, character(1)))
}

# The number of samples with a detected peak, and the median m/z and rt of the
# peaks when the chromatograms carry none.
.add_peak_summaries <- function(index, peaks) {
    if (nrow(peaks) == 0 || !"feature_id" %in% colnames(peaks)) {
        index$n_samples <- rep(0L, nrow(index))
        return(index)
    }

    rows <- lapply(index$feature_ids, function(ids) {
        peaks[peaks$feature_id %in% ids, , drop = FALSE]
    })
    median_of <- function(column) {
        vapply(rows, function(df) {
            if (nrow(df) == 0) NA_real_ else stats::median(df[[column]])
        }, numeric(1))
    }
    if (!"mz" %in% colnames(index)) index$mz <- median_of("mz")
    if (!"rt" %in% colnames(index)) index$rt <- median_of("rt")
    index$n_samples <- vapply(rows, function(df) {
        length(unique(df$sample_index))
    }, integer(1))
    index
}
