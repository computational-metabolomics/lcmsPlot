#' Save a figure with its provenance
#'
#' `lp_save()` saves the plot of a pipeline with [ggplot2::ggsave()] and
#' records how it was made next to it, in a `.provenance.json` file named
#' after it (`figure.pdf` gives `figure.pdf.provenance.json`): the versions of
#' lcmsPlot, R, Bioconductor and every loaded package, the layers of the
#' pipeline with their arguments, and the input files with their size and
#' modification time. With `data = TRUE`, the data drawn in the figure is also
#' written by [lp_export_data()] to a folder of the same name
#' (`figure_data/`).
#'
#' @param object An instance of class `lcmsPlotClass`, with data to plot.
#' @param file The path of the figure. Its extension sets the format (`.pdf`,
#' `.png`, `.svg`, `.tiff`, ...), as in [ggplot2::ggsave()].
#' @param width,height The size of the figure, in `units`.
#' @param units The units of `width` and `height`: `"in"`, `"cm"`, `"mm"`, or
#' `"px"`.
#' @param dpi The resolution of raster formats.
#' @param data A `logical` value indicating whether to also write the plotted
#' data (see [lp_export_data()]).
#' @param provenance A `logical` value indicating whether to write the
#' provenance file.
#' @param checksums A `logical` value indicating whether to record the MD5
#' checksum of every input file in the provenance. Reading large raw files can
#' take a while, so it is off by default.
#' @param ... Other arguments passed to [ggplot2::ggsave()].
#' @return The path of the figure (invisibly).
#' @export
#' @examples
#' raw_files <- dir(
#'    system.file("cdf", package = "faahKO"),
#'    full.names = TRUE,
#'    recursive = TRUE)[1:2]
#'
#' p <- lcmsPlot(raw_files) +
#'     lp_chromatogram(features = rbind(c(
#'         mzmin = 334.9, mzmax = 335.1, rtmin = 2700, rtmax = 2900))) +
#'     lp_arrange(group_by = "sample_id")
#'
#' file <- tempfile(fileext = ".pdf")
#' lp_save(p, file, data = TRUE)
#'
#' ## The provenance and the plotted data, next to the figure
#' prov <- jsonlite::read_json(paste0(file, ".provenance.json"))
#' prov$pipeline$layers[[1]]$layer
#' dir(sub("\\.pdf$", "_data", file))
lp_save <- function(
    object,
    file,
    width = 7,
    height = 5,
    units = "in",
    dpi = 300,
    data = FALSE,
    provenance = TRUE,
    checksums = FALSE,
    ...
) {
    .check_plot_data(object, "lp_save")
    object <- .render_plot(object, .additional_datasets(object))

    format <- tolower(tools::file_ext(file))
    args <- list(
        filename = file, plot = object@plot, width = width,
        height = height, units = units, dpi = dpi, ...)
    if (format == "pdf" && is.null(args$title)) {
        args$title <- "lcmsPlot figure"
    }
    do.call(ggsave, args)

    data_dir <- NULL
    if (isTRUE(data)) {
        data_dir <- paste0(tools::file_path_sans_ext(file), "_data")
        lp_export_data(object, data_dir)
    }
    if (isTRUE(provenance)) {
        output <- list(
            file = basename(file), format = format, width = width,
            height = height, units = units, dpi = dpi,
            data = if (!is.null(data_dir)) basename(data_dir))
        .write_provenance(.provenance(object, checksums, output), file)
    }

    invisible(file)
}

#' Export the plotted data in an open format
#'
#' `lp_export_data()` writes the data a pipeline draws as plain tables, so
#' that a figure can be checked, re-plotted, or reused without R: one file per
#' dataset of the plot (`chromatograms`, `spectra`, `mass_traces`, ...), in
#' the tidy form the plot is drawn from, with one row per data point and the
#' sample metadata joined in. The folder also holds:
#'
#' - `detected_peaks`: the detected peaks, when there are any;
#' - `samples`: the sample metadata;
#' - `columns`: a data dictionary of every table, with the type, unit and
#'   description of each column.
#'
#' Retention times are in seconds, whatever the unit of the plot axes.
#'
#' @param object An instance of class `lcmsPlotClass`, with data to plot.
#' @param dir The folder to write the tables to. It is created if needed.
#' @param format The format of the tables: `"csv"` (comma-separated) or
#' `"tsv"` (tab-separated).
#' @return A named `character` vector of the paths of the written files
#' (invisibly).
#' @export
#' @examples
#' raw_files <- dir(
#'    system.file("cdf", package = "faahKO"),
#'    full.names = TRUE,
#'    recursive = TRUE)[1:2]
#'
#' p <- lcmsPlot(raw_files) +
#'     lp_chromatogram(features = rbind(c(
#'         mzmin = 334.9, mzmax = 335.1, rtmin = 2700, rtmax = 2900)))
#'
#' paths <- lp_export_data(p, file.path(tempdir(), "plot-data"))
#' basename(paths)
#' head(read.csv(paths[["chromatograms"]]))
lp_export_data <- function(object, dir, format = c("csv", "tsv")) {
    .check_plot_data(object, "lp_export_data")
    format <- match.arg(format)

    tables <- .plot_datasets(object, .additional_datasets(object))
    peaks <- .supporting_datasets(object)$detected_peaks
    if (nrow(peaks) > 0) {
        tables$detected_peaks <- peaks
    }
    tables$samples <- sampleMetadata(object)
    tables$columns <- .column_dictionary(tables)

    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    paths <- vapply(names(tables), function(name) {
        path <- file.path(dir, paste0(name, ".", format))
        data.table::fwrite(
            tables[[name]], path, sep = if (format == "tsv") "\t" else ",",
            na = "NA")
        path
    }, character(1))

    invisible(paths)
}

.check_plot_data <- function(object, fn) {
    if (!is(object, "lcmsPlotClass") || !.has_data(object)) {
        stop(
            fn, ": 'object' must be an lcmsPlot object with data; ",
            "add a layer such as lp_chromatogram() first.", call. = FALSE)
    }
}

# The units and descriptions of the columns lcmsPlot writes. Columns not
# listed here (such as sample metadata) are described as NA.
.known_columns <- tibble::tribble(
    ~column, ~unit, ~description,
    "rt", "s", "Retention time",
    "rtmin", "s", "Lower bound of the retention time window",
    "rtmax", "s", "Upper bound of the retention time window",
    "rt_raw", "s", "Retention time before alignment",
    "rt_adj", "s", "Retention time after alignment",
    "diff", "s", "Adjusted minus raw retention time",
    "mz", "m/z", "Mass-to-charge ratio",
    "mzmin", "m/z", "Lower bound of the m/z window",
    "mzmax", "m/z", "Upper bound of the m/z window",
    "precursor_mz", "m/z", "Precursor m/z of an MS2 spectrum",
    "isolation_window_half_width", "m/z",
    "Half width of the precursor isolation window",
    "intensity", "a.u.", "Signal intensity",
    "into", "a.u.", "Integrated intensity of the detected peak",
    "maxo", "a.u.", "Maximum intensity of the detected peak",
    "sample_index", NA, "Row number of the sample in the samples table",
    "sample_id", NA, "Sample identifier",
    "sample_path", NA, "Path of the sample's data file",
    "feature_id", NA, "Feature identifier",
    "name", NA, "Compound name",
    "reference", NA, "Whether the spectrum is a library reference",
    "in_purity", NA, "Whether the spectrum is used for precursor purity"
)

# The data dictionary of the exported tables: one row per column of each.
.column_dictionary <- function(tables) {
    rows <- lapply(names(tables), function(name) {
        df <- tables[[name]]
        tibble(
            table = name,
            column = colnames(df),
            type = vapply(df, function(col) class(col)[1], character(1)))
    })
    left_join(
        bind_rows(rows), .known_columns, by = "column")
}
