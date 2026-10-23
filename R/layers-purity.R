#' Overlay precursor ion purity scores on a chromatogram
#'
#' `lp_purity_overlay()` adds a `geom_point` layer to the chromatogram panel
#' where each triangle marks the retention time of an MS/MS fragmentation event,
#' coloured by the interpolated precursor ion purity (`inPurity`). Requires the
#' data object to be a `purityA` result from the msPurity package and
#' `lp_chromatogram()` to be called first.
#'
#' @param sample_ids A `character` vector of sample IDs to include.
#' `NULL` uses all samples.
#' @param threshold A `numeric` value in `[0, 1]` drawn as a reference
#' on the colour scale midpoint. `NULL` defaults to `0.5`.
#' @param point_size A `numeric` value controlling the point size.
#' @return A layer function for use with the `+` operator.
#' @export
#' @examplesIf requireNamespace("msPurity", quietly = TRUE)
#' ## lcmsPlot ships two DDA files carrying real precursor isolation metadata.
#' mzml_dir <- file.path(tempdir(), "purity-example")
#' utils::unzip(
#'     system.file("extdata", "standards-mzml.zip", package = "lcmsPlot"),
#'     exdir = mzml_dir)
#'
#' ms2_files <- file.path(
#'     mzml_dir, "mzml",
#'     c("l-proline-MS2.mzML", "l-kynurenine-MS2.mzML"))
#'
#' pa <- msPurity::purityA(ms2_files)
#'
#' ## The overlay annotates an existing chromatogram, so the retention time
#' ## window has to span the fragmentation events to show anything.
#' features <- data.frame(
#'     sample_id = c("l-proline-MS2", "l-kynurenine-MS2"),
#'     mz = c(116.0703793, 209.091824),
#'     rt = c(425.74, 365.72))
#'
#' p <- lcmsPlot(pa) +
#'     lp_chromatogram(features = features, ppm = 10, rt_tol = 60) +
#'     lp_purity_overlay(threshold = 0.7)
#' p
lp_purity_overlay <- function(
        sample_ids = NULL,
        threshold = NULL,
        point_size = 2
) {
    args <- as.list(environment())

    make_interface_function(
        name = "lp_purity_overlay",
        args_list = args,
        fn = function(obj) {
            .require_source(
                obj, is_mspurity_data, "lp_purity_overlay requires a purityA data object")

            if (!isTRUE(obj@options$chromatograms$show)) {
                stop("lp_purity_overlay must be called after lp_chromatogram()")
            }

            obj@options$purity_overlay <- .layer_options(args, show = TRUE)
            obj@options$purity_scores_sample_ids <- sample_ids

            obj@data <- create_purity_scores(obj@data, obj@options)
            return(obj)
        }
    )
}

#' Visualise the isolation window for a precursor fragmentation event
#'
#' `lp_isolation_window()` extracts the MS1 survey spectrum at the scan
#' immediately preceding a chosen MS2 event from a `purityA` object and
#' annotates it with a semi-transparent rectangle spanning the isolation window
#' and a dashed line at the precursor m/z. The `inPurity` score is shown as a
#' plot subtitle.
#'
#' @param pid An integer or integer vector of precursor IDs (the `pid` column in
#' `purityA@@puritydf`) selecting which event(s) to display. `NULL` shows all
#' events (use with care on large datasets).
#' @param sample_id A `character` sample ID to filter by when `pid` is `NULL`.
#' @param half_width A `numeric` value (in Da) for the isolation window
#' half-width on each side of the precursor m/z. Default is `0.5`.
#' @param zoom_factor A `numeric` multiplier controlling how far beyond the
#' isolation window the x-axis extends. The visible range is
#' `precursor_mz ± zoom_factor * half_width`. Default is `3`.
#' @return A layer function for use with the `+` operator.
#' @export
#' @examplesIf requireNamespace("msPurity", quietly = TRUE)
#' mzml_dir <- file.path(tempdir(), "purity-example")
#' utils::unzip(
#'     system.file("extdata", "standards-mzml.zip", package = "lcmsPlot"),
#'     exdir = mzml_dir)
#'
#' ms2_files <- file.path(
#'     mzml_dir, "mzml",
#'     c("l-proline-MS2.mzML", "l-kynurenine-MS2.mzML"))
#'
#' pa <- msPurity::purityA(ms2_files)
#'
#' ## Inspect the fragmentation events available to plot.
#' slot(pa, "puritydf")[, c("pid", "precursorMZ", "precursorRT", "inPurity")]
#'
#' ## `pid = 4` is a co-isolated precursor (inPurity ~ 0.51). These files
#' ## record an isolation width of 0.6 Da on each side of the precursor.
#' p <- lcmsPlot(pa) +
#'     lp_isolation_window(pid = 4, half_width = 0.6)
#' p
lp_isolation_window <- function(
        pid = NULL,
        sample_id = NULL,
        half_width = 0.5,
        zoom_factor = 3
) {
    args <- as.list(environment())

    make_interface_function(
        name = "lp_isolation_window",
        args_list = args,
        fn = function(obj) {
            .require_source(
                obj, is_mspurity_data, "lp_isolation_window requires a purityA data object")

            obj@options$isolation_window <- .layer_options(args, show = TRUE)
            obj@options$spectra$show <- TRUE

            obj@data <- create_isolation_window(
                obj@data, obj@options, pid, sample_id)
            return(obj)
        }
    )
}

#' Plot precursor ion purity scores as a timeline
#'
#' `lp_purity_timeline()` generates a scatter plot of `inPurity` (y-axis)
#' versus retention time (x-axis), one point per MS/MS acquisition event,
#' coloured by sample. An optional horizontal dashed line marks a purity
#' threshold. Requires the data object to be a `purityA` result.
#'
#' @param sample_ids A `character` vector of sample IDs to include.
#' `NULL` uses all samples.
#' @param threshold A `numeric` value in `[0, 1]` drawn as a horizontal
#' reference line. `NULL` suppresses the line.
#' @return A layer function for use with the `+` operator.
#' @export
#' @examplesIf requireNamespace("msPurity", quietly = TRUE)
#' mzml_dir <- file.path(tempdir(), "purity-example")
#' utils::unzip(
#'     system.file("extdata", "standards-mzml.zip", package = "lcmsPlot"),
#'     exdir = mzml_dir)
#'
#' ms2_files <- file.path(
#'     mzml_dir, "mzml",
#'     c("l-proline-MS2.mzML", "l-kynurenine-MS2.mzML"))
#'
#' pa <- msPurity::purityA(ms2_files)
#'
#' ## One point per MS/MS event, coloured by sample.
#' p <- lcmsPlot(pa) +
#'     lp_purity_timeline(threshold = 0.7)
#' p
lp_purity_timeline <- function(
        sample_ids = NULL,
        threshold = NULL
) {
    args <- as.list(environment())

    make_interface_function(
        name = "lp_purity_timeline",
        args_list = args,
        fn = function(obj) {
            .require_source(
                obj, is_mspurity_data, "lp_purity_timeline requires a purityA data object")

            obj@options$purity_timeline <- .layer_options(args, show = TRUE)
            obj@options$purity_scores_sample_ids <- sample_ids

            obj@data <- create_purity_scores(obj@data, obj@options)
            return(obj)
        }
    )
}

#' Plot the distribution of precursor ion purity scores per sample
#'
#' `lp_purity_distribution()` generates a violin (or box/jitter) plot of
#' `inPurity` scores grouped by sample. An optional horizontal dashed line
#' marks a purity acceptance threshold. Requires the data object to be a
#' `purityA` result from the msPurity package.
#'
#' @param sample_ids A `character` vector of sample IDs to include.
#' `NULL` uses all samples.
#' @param threshold A `numeric` value in `[0, 1]` drawn as a horizontal
#' reference line. `NULL` suppresses the line.
#' @param type A `character` value; one of `"violin"`, `"boxplot"`,
#' or `"jitter"`. Defaults to `"violin"`.
#' @return A layer function for use with the `+` operator.
#' @export
#' @examplesIf requireNamespace("msPurity", quietly = TRUE)
#' mzml_dir <- file.path(tempdir(), "purity-example")
#' utils::unzip(
#'     system.file("extdata", "standards-mzml.zip", package = "lcmsPlot"),
#'     exdir = mzml_dir)
#'
#' ms2_files <- file.path(
#'     mzml_dir, "mzml",
#'     c("l-proline-MS2.mzML", "l-kynurenine-MS2.mzML"))
#'
#' pa <- msPurity::purityA(ms2_files)
#'
#' ## Compare the spread of purity scores between samples.
#' p <- lcmsPlot(pa) +
#'     lp_purity_distribution(threshold = 0.7, type = "boxplot")
#' p
lp_purity_distribution <- function(
        sample_ids = NULL,
        threshold = NULL,
        type = "violin"
) {
    args <- as.list(environment())

    make_interface_function(
        name = "lp_purity_distribution",
        args_list = args,
        fn = function(obj) {
            .require_source(
                obj, is_mspurity_data, "lp_purity_distribution requires a purityA data object")

            obj@options$purity_distribution <- .layer_options(args, show = TRUE)
            obj@options$purity_scores_sample_ids <- sample_ids

            obj@data <- create_purity_scores(obj@data, obj@options)
            return(obj)
        }
    )
}
