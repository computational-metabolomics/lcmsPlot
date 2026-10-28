# The datasets drawn as panels, in panel order.
DATASET_TYPES <- c(
    "chromatograms",
    "mass_traces",
    "spectra",
    "peak_density",
    "peak_count_image",
    "intensity_maps",
    "total_ion_current",
    "rt_diff"
)

#' Convert a dataset into the tidy form used for plotting
#'
#' `as_plot_df()` is the single place where the data structures of an
#' `lcmsPlotDataContainer` are turned into the long `tibble`s the plotting
#' functions draw from. The conversion reads every chromatogram's or
#' spectrum's peaks in one call and builds the columns as whole vectors.
#'
#' - `Chromatograms`: one row per data point with `rt` and `intensity`, plus
#'   the annotation columns of `chromData()` (such as `sample_index`,
#'   `feature_id`, or `name`). Window bounds and empty core variables are
#'   left out.
#' - `Spectra`: one row per peak with `mz`, `intensity`, and `rt`, plus the
#'   spectra variables set by lcmsPlot (such as `sample_index`, `feature_id`,
#'   and `reference`).
#' - `data.frame`: returned as a `tibble`.
#'
#' @param x A `Chromatograms`, `Spectra`, or `data.frame` object.
#' @param ... Ignored.
#' @return A `tibble`.
#' @keywords internal
setGeneric("as_plot_df", function(x, ...) standardGeneric("as_plot_df"))

#' @rdname as_plot_df
setMethod("as_plot_df", "Chromatograms", function(x, ...) {
    if (length(x) == 0) {
        return(tibble())
    }

    peaks <- Chromatograms::peaksData(x, columns = c("rtime", "intensity"))
    n_points <- vapply(peaks, nrow, integer(1))
    chrom_data <- as.data.frame(Chromatograms::chromData(x))
    annotations <- .chrom_annotation_columns(chrom_data)

    bind_cols(
        tibble(
            rt = unlist(lapply(peaks, `[[`, "rtime"), use.names = FALSE),
            intensity = unlist(
                lapply(peaks, `[[`, "intensity"), use.names = FALSE)
        ),
        as_tibble(chrom_data[
            rep.int(seq_len(nrow(chrom_data)), n_points),
            annotations,
            drop = FALSE])
    )
})

# The chromData() columns worth plotting: lcmsPlot's annotations, without the
# extraction window bounds and the core variables, except an `mz` that holds a
# value.
.chrom_annotation_columns <- function(chrom_data) {
    technical <- c(
        names(Chromatograms::coreChromVariables()),
        "rtMin", "rtMax", "chromSpectraIndex")
    keep <- setdiff(colnames(chrom_data), technical)
    if ("mz" %in% colnames(chrom_data) && !all(is.na(chrom_data$mz))) {
        keep <- c(keep, "mz")
    }
    keep
}

# Spectra variables lcmsPlot sets on the spectra it extracts.
.plot_spectra_variables <- c(
    "sample_index",
    "feature_id",
    "reference",
    "in_purity",
    "precursor_mz",
    "isolation_window_half_width"
)

#' @rdname as_plot_df
setMethod("as_plot_df", "Spectra", function(x, ...) {
    if (length(x) == 0) {
        return(tibble())
    }

    peaks <- Spectra::peaksData(x, return.type = "list")
    n_peaks <- vapply(peaks, nrow, integer(1))
    points <- do.call(rbind, peaks)
    variables <- intersect(
        .plot_spectra_variables, Spectra::spectraVariables(x))
    spectra_data <- as.data.frame(
        Spectra::spectraData(x, c("rtime", variables)))
    rows <- rep.int(seq_len(nrow(spectra_data)), n_peaks)

    bind_cols(
        tibble(
            mz = as.numeric(points[, 1]),
            intensity = as.numeric(points[, 2]),
            rt = spectra_data$rtime[rows]
        ),
        as_tibble(spectra_data[rows, variables, drop = FALSE])
    )
})

#' @rdname as_plot_df
setMethod("as_plot_df", "data.frame", function(x, ...) {
    as_tibble(x)
})

# The tidy form of one dataset of a container.
.plot_dataset <- function(container, name) {
    switch(
        name,
        chromatograms = as_plot_df(container@chromatograms),
        spectra = as_plot_df(container@spectra),
        as_plot_df(.summary(container, name))
    )
}

.render_plot <- function(object, additional_datasets) {
    object@plot <- plot_data(
        .plot_datasets(object, additional_datasets), object)

    return(object)
}

# The datasets drawn as panels, in panel order, with the sample metadata
# joined in. These are what the plot shows and what lp_export_data() writes.
.plot_datasets <- function(object, additional_datasets) {
    dataset_types <- c(DATASET_TYPES, names(additional_datasets))
    sample_table <- .sample_table(object@data)

    datasets <- lapply(dataset_types, function(dataset_name) {
        if (!isTRUE(object@options[[dataset_name]]$show)) {
            return(NULL)
        }

        data_df <- if (dataset_name %in% names(additional_datasets)) {
            additional_datasets[[dataset_name]]
        } else {
            .plot_dataset(object@data, dataset_name)
        }

        if (nrow(data_df) == 0) {
            warning("Empty dataset ", dataset_name)
            return(NULL)
        }

        .join_sample_data(data_df, sample_table)
    })
    names(datasets) <- dataset_types
    remove_null_elements(datasets)
}

.has_data <- function(object) {
    data <- object@data
    has_summary <- any(vapply(data@summaries, nrow, integer(1)) > 0)

    # A standalone peak-boundary panel is drawn from the detected peaks, which
    # are not a dataset of their own. Detected peaks alone are not a plot
    # request, so the option has to be on as well.
    length(data@chromatograms) > 0 ||
        length(data@spectra) > 0 ||
        has_summary ||
        (isTRUE(object@options$chrom_peak_rects$show) &&
            nrow(data@detected_peaks) > 0)
}

# Datasets that are derived at render time rather than stored on their own.
# `.render_plot()` appends these to DATASET_TYPES, so each name needs a
# matching options entry carrying a logical `show`.
.additional_datasets <- function(object) {
    additional <- list()
    purity_scores <- .summary(object@data, "purity_scores")

    if (nrow(purity_scores) > 0) {
        if (isTRUE(object@options$purity_timeline$show)) {
            additional$purity_timeline <- purity_scores
        }
        if (isTRUE(object@options$purity_distribution$show)) {
            additional$purity_distribution <- purity_scores
        }
    }

    # The peak rectangles own a panel only when there is no host layer to
    # decorate. Deriving that here rather than at `+` time keeps it independent
    # of the order the layers were added in.
    has_host <- isTRUE(object@options$intensity_maps$show) ||
        isTRUE(object@options$mass_traces$show)
    peaks <- object@data@detected_peaks

    if (isTRUE(object@options$chrom_peak_rects$show) && !has_host &&
        nrow(peaks) > 0) {
        additional$chrom_peak_rects <- peaks |>
            select(all_of(c(
                "mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax",
                "sample_index")))
    }

    additional
}

# The datasets that decorate the panels rather than being drawn as panels of
# their own, with the sample metadata joined in.
.supporting_datasets <- function(object) {
    sample_table <- .sample_table(object@data)
    list(
        detected_peaks = .join_sample_data(
            object@data@detected_peaks, sample_table),
        purity_scores = .summary(object@data, "purity_scores")
    )
}
