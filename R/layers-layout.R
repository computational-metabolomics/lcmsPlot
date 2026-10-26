#' Define the arrangement of chromatograms
#'
#' The `lp_arrange` function specifies how chromatograms should be arranged
#' when visualised. It determines the grouping metadata factor through
#' the `group_by` parameter.
#'
#' @param group_by A `character` value determining the
#' column to group by in the samples metadata.
#' @return A function that takes an `lcmsPlot` object and returns a modified
#' version with the specified arrangement options
#' stored in `options$arrangement`. It is intended for use with the `+`
#' operator, which incrementally layers new data or visual components
#' onto the `lcmsPlot` object.
#' @export
#' @examples
#' raw_files <- dir(
#'    system.file("cdf", package = "faahKO"),
#'    full.names = TRUE,
#'    recursive = TRUE)[1:5]
#'
#' ## Plots chromatograms overlayed without specifying a grouping factor
#' p <- lcmsPlot(raw_files) +
#'   lp_chromatogram(aggregation_fun = "max")
#' p
#'
#' ## Plots chromatograms overlayed specifying a grouping factor
#' ## (e.g., sample_id)
#' p <- p + lp_arrange(group_by = "sample_id")
#' p
lp_arrange <- function(group_by) {
    args <- as.list(environment())

    make_interface_function(
        name = "lp_arrange",
        args_list = args,
        fn = function(obj) {
            obj@options$arrangement <- args
            return(obj)
        }
    )
}

#' Define the plot's faceting
#'
#' The `lp_facets` function arranges plots into a grid based
#' on a metadata factor, creating a series of smaller plots (facets).
#'
#' @param facets A `character` vector of factors from the sample metadata to use
#' for faceting.
#' @param ncol A `numeric` value indicating the number of columns in the layout.
#' @param nrow A `numeric` value indicating the number of rows in the layout.
#' @param free_x A `logical` value indicating whether the x-axis scales
#' are allowed to vary across panels.
#' @param free_y A `logical` value indicating whether the y-axis scales
#' are allowed to vary across panels.
#' @return A function that takes an `lcmsPlot` object and returns a modified
#' version with the specified faceting options stored in `options$facets`.
#' It is intended for use with the `+` operator, which incrementally layers
#' new data or visual components onto the `lcmsPlot` object.
#' @export
#' @examples
#' raw_files <- dir(
#'    system.file("cdf", package = "faahKO"),
#'    full.names = TRUE,
#'    recursive = TRUE)[1:5]
#'
#' ## Plots chromatograms overlayed
#' p <- lcmsPlot(raw_files) +
#'   lp_chromatogram(aggregation_fun = "max")
#' p
#'
#' ## Using lp_facets we create facets for each sample_id
#' p <- p + lp_facets(facets = "sample_id")
#' p
lp_facets <- function(
    facets,
    ncol = NULL,
    nrow = NULL,
    free_x = FALSE,
    free_y = FALSE
) {
    args <- as.list(environment())

    make_interface_function(
        name = "lp_facets",
        args_list = args,
        fn = function(obj) {
            obj@options$facets <- args
            return(obj)
        }
    )
}

#' Define a gridded plot
#'
#' The `lp_grid` function arranges plots into a matrix of panels defined
#' by row and column faceting metadata factors.
#'
#' @param rows A `character` value indicating the factors that
#' represent rows.
#' @param cols A `character` value indicating the factors that
#' represent columns.
#' @param free_x A `logical` value indicating whether the x-axis scales
#' are allowed to vary across panels.
#' @param free_y A `logical` value indicating whether the y-axis scales
#' are allowed to vary across panels.
#' @return A function that takes an `lcmsPlot` object and returns a modified
#' version with the specified grid options stored in `options$grid`.
#' It is intended for use with the `+` operator, which incrementally layers
#' new data or visual components onto the `lcmsPlot` object.
#' @export
#' @examples
#' raw_files <- dir(
#'   system.file("cdf", package = "faahKO"),
#'   full.names = TRUE,
#'   recursive = TRUE
#' )[1:4]
#'
#' ## Create metadata for the samples
#' metadata <- data.frame(
#'   sample_id = sub("\\.CDF", "", basename(raw_files)),
#'   factor1 = c("S", "S", "C", "C"),
#'   factor2 = c("T", "U", "T", "U")
#' )
#'
#' ## Create feature chromatograms for the specified samples
#' p <- lcmsPlot(raw_files, metadata = metadata) +
#'   lp_chromatogram(features = rbind(c(
#'     mzmin = 334.9,
#'     mzmax = 335.1,
#'     rtmin = 2700,
#'     rtmax = 2900)))
#' p
#'
#' ## Arrange chromatograms in a grid split by experimental factors
#' ## Rows correspond to `factor1` and columns correspond to `factor2`
#' p <- p + lp_grid(rows = "factor1", cols = "factor2")
#' p
lp_grid <- function(rows, cols, free_x = FALSE, free_y = FALSE) {
    args <- as.list(environment())

    make_interface_function(
        name = "lp_grid",
        args_list = args,
        fn = function(obj) {
            obj@options$grid <- args
            return(obj)
        }
    )
}

#' Define the labels of the plot
#'
#' The `lp_labels` function allows the specification
#' of the plot title and the legend title.
#'
#' @param title A `character` value indicating the plot title.
#' @param legend A `character` value indicating the legend's title.
#' @return A function that takes an `lcmsPlot` object and returns a modified
#' version with the specified label options stored in `options$labels`.
#' It is intended for use with the `+` operator, which incrementally layers
#' new data or visual components onto the `lcmsPlot` object.
#' @export
#' @examples
#' raw_files <- dir(
#'    system.file("cdf", package = "faahKO"),
#'    full.names = TRUE,
#'    recursive = TRUE)[1:5]
#'
#' ## Create a chromatogram plot by grouping samples into batches
#' ## By default, the legend is derived from the grouping variable
#' p <- lcmsPlot(raw_files, batch_size = 2) +
#'   lp_chromatogram(features = rbind(c(
#'     mzmin = 334.9,
#'     mzmax = 335.1,
#'     rtmin = 2700,
#'     rtmax = 2900))) +
#'   lp_arrange(group_by = "sample_id")
#' p
#'
#' ## Customise the legend label
#' p <- p + lp_labels(legend = "Sample")
#' p
lp_labels <- function(title = NULL, legend = NULL) {
    args <- as.list(environment())

    make_interface_function(
        name = "lp_labels",
        args_list = args,
        fn = function(obj) {
            obj@options$labels <- args
            return(obj)
        }
    )
}

#' Define the legend layout
#'
#' @param position A `character` value indicating the legend's position.
#' One of `"top"`, `"right"`, `"bottom"`, `"left"`, or `"inside"`.
#' @return A function that takes an `lcmsPlot` object and returns a modified
#' version with the specified legend options stored in `options$legend`.
#' It is intended for use with the `+` operator, which incrementally layers
#' new data or visual components onto the `lcmsPlot` object.
#' @export
#' @examples
#' raw_files <- dir(
#'    system.file("cdf", package = "faahKO"),
#'    full.names = TRUE,
#'    recursive = TRUE)[1:5]
#'
#' ## Create a chromatogram plot grouped by sample with a custom legend label
#' p <- lcmsPlot(raw_files, batch_size = 2) +
#'   lp_chromatogram(features = rbind(c(
#'     mzmin = 334.9,
#'     mzmax = 335.1,
#'     rtmin = 2700,
#'     rtmax = 2900))) +
#'   lp_arrange(group_by = "sample_id") +
#'   lp_labels(legend = "Sample")
#' p
#'
#' ## Move the legend below the plot
#' p <- p + lp_legend(position = "bottom")
#' p
lp_legend <- function(position = NULL)  {
    args <- as.list(environment())

    make_interface_function(
        name = "lp_legend",
        args_list = args,
        fn = function(obj) {
            obj@options$legend <- args
            return(obj)
        }
    )
}

#' Define the theme and colour palette of the plot
#'
#' The `lp_theme` function styles every panel of the plot with a theme preset
#' and colours the series (the `group_by` column of [lp_arrange()], and the
#' highlighted peaks) with a discrete palette. Continuous colour scales, such
#' as those of intensity maps or purity scores, are left unchanged.
#'
#' A theme added after [lp_get_plot()] (e.g. `lp_get_plot() + theme(...)`)
#' still overrides the preset.
#'
#' @param preset The theme preset. One of `"minimal"` (default),
#' `"classic"`, `"publication"` (white background, axis lines, no minor grid,
#' and an 8 pt base size for print), or any `ggplot2` `theme` object.
#' @param palette The discrete palette. One of `"okabe_ito"` (default;
#' colourblind-safe, falling back to viridis beyond 8 groups), `"viridis"`,
#' `"ggplot2"` (the default ggplot2 hues), or a `character` vector of colours.
#' A named vector maps each group to its colour; an unnamed one is used in
#' order.
#' @param base_size A `numeric` value indicating the base font size in points.
#' `NULL` uses the preset's default (11, or 8 for `"publication"`).
#' @param line_width A `numeric` value indicating the width of the plotted
#' lines (chromatogram traces, spectrum sticks) and of the theme's lines.
#' `NULL` scales it with `base_size`. `base_size` and `line_width` only affect
#' the theme when `preset` is one of the named presets.
#' @return A function that takes an `lcmsPlot` object and returns a modified
#' version with the specified theme options stored in `options$theme`.
#' It is intended for use with the `+` operator, which incrementally layers
#' new data or visual components onto the `lcmsPlot` object.
#' @export
#' @examples
#' raw_files <- dir(
#'    system.file("cdf", package = "faahKO"),
#'    full.names = TRUE,
#'    recursive = TRUE)[1:4]
#'
#' p <- lcmsPlot(raw_files) +
#'   lp_chromatogram(features = rbind(c(
#'     mzmin = 334.9,
#'     mzmax = 335.1,
#'     rtmin = 2700,
#'     rtmax = 2900))) +
#'   lp_arrange(group_by = "sample_id")
#'
#' ## Style the plot for print
#' p + lp_theme(preset = "publication", palette = "okabe_ito", base_size = 8)
#'
#' ## Colour each sample explicitly
#' p + lp_theme(palette = c(
#'   ko15 = "#0072B2", ko16 = "#56B4E9",
#'   ko18 = "#D55E00", ko19 = "#E69F00"))
lp_theme <- function(
    preset = "minimal",
    palette = "okabe_ito",
    base_size = NULL,
    line_width = NULL
) {
    .check_theme_preset(preset)
    .check_palette(palette)
    .check_positive_number(base_size, "base_size")
    .check_positive_number(line_width, "line_width")

    args <- as.list(environment())

    make_interface_function(
        name = "lp_theme",
        args_list = args,
        fn = function(obj) {
            obj@options$theme <- args
            return(obj)
        }
    )
}

.check_theme_preset <- function(preset) {
    if (inherits(preset, "theme")) {
        return(invisible())
    }
    if (!is.character(preset) || length(preset) != 1 ||
        !preset %in% .theme_presets) {
        stop(
            "lp_theme: 'preset' must be one of ",
            paste0('"', .theme_presets, '"', collapse = ", "),
            ", or a ggplot2 theme.", call. = FALSE)
    }
}

.check_palette <- function(palette) {
    if (!is.character(palette) || length(palette) == 0) {
        stop("lp_theme: 'palette' must be a character vector.", call. = FALSE)
    }
    if (length(palette) == 1 && is.null(names(palette)) &&
        palette %in% .palette_presets) {
        return(invisible())
    }

    valid <- vapply(palette, function(colour) {
        tryCatch({
            grDevices::col2rgb(colour)
            TRUE
        }, error = function(e) FALSE)
    }, logical(1))
    if (!all(valid)) {
        stop(
            "lp_theme: 'palette' must be one of ",
            paste0('"', .palette_presets, '"', collapse = ", "),
            ", or a vector of colours. Not colours: ",
            paste0('"', palette[!valid], '"', collapse = ", "),
            call. = FALSE)
    }
}

.check_positive_number <- function(x, name) {
    if (is.null(x)) {
        return(invisible())
    }
    if (!is.numeric(x) || length(x) != 1 || is.na(x) || x <= 0) {
        stop(
            "lp_theme: '", name, "' must be a single positive number.",
            call. = FALSE)
    }
}

#' Define a vertical line on a retention time value
#'
#' @param intercept A `numeric` value indicating the retention time
#' axis (x-axis) intercept.
#' @param line_type A `character` value indicating the line type.
#' One of `"solid"`, `"dashed"`, `"dotted"`,
#' `"dotdash"`, `"longdash"`, `"twodash"`.
#' @param color A `character` value indicating the line color.
#' @return A function that takes an `lcmsPlot` object and returns a modified
#' version with the specified RT line options stored in `options$rt_lines`.
#' It is intended for use with the `+` operator, which incrementally layers
#' new data or visual components onto the `lcmsPlot` object.
#' @export
#' @examples
#' raw_files <- dir(
#'    system.file("cdf", package = "faahKO"),
#'    full.names = TRUE,
#'    recursive = TRUE)[1:4]
#'
#' ## Create chromatogram plots faceted by sample
#' p <- lcmsPlot(raw_files) +
#'   lp_chromatogram(features = rbind(c(
#'     mzmin = 334.9,
#'     mzmax = 335.1,
#'     rtmin = 2700,
#'     rtmax = 2900))) +
#'   lp_facets(facets = 'sample_id', ncol = 4)
#' p
#'
#' ## Add a vertical retention time reference line
#' p <- p + lp_rt_line(intercept = 2800, line_type = 'solid', color = 'red')
#' p
lp_rt_line <- function(intercept, line_type = 'dashed', color = 'black') {
    args <- as.list(environment())

    make_interface_function(
        name = "lp_rt_line",
        args_list = args,
        fn = function(obj) {
            obj@options$rt_lines <- append(obj@options$rt_lines, list(args))
            return(obj)
        }
    )
}

#' Define the plot layout
#'
#' @param design Specification of the location of areas in the layout
#' See https://patchwork.data-imaginist.com/reference/wrap_plots.html
#' @return A function that takes an `lcmsPlot` object and returns a modified
#' version with the specified layout options stored in `options$layout`.
#' It is intended for use with the `+` operator, which incrementally layers
#' new data or visual components onto the `lcmsPlot` object.
#' @export
#' @examples
#' data_obj <- get_XCMSnExp_object_example(indices = 1)
#'
#' ## Plot chromatograms and spectra for selected samples and features
#' p <- lcmsPlot(data_obj, sample_id_column = 'sample_name') +
#'   lp_chromatogram(
#'     features = rbind(
#'       c(mzmin = 334.9, mzmax = 335.1, rtmin = 2700, rtmax = 2900),
#'       c(mzmin = 278.99721, mzmax = 279.00279, rtmin = 2740, rtmax = 2840)
#'     ),
#'     sample_ids = 'ko15',
#'     highlight_peaks = TRUE
#'   ) +
#'   lp_spectra(mode = "closest_apex", ms_level = 1) +
#'   lp_facets(facets = "feature_id", ncol = 2)
#'
#' ## Customise panel layout to place chromatogram above spectra
#' p <- p + lp_layout(design = "C\nS\nS")
lp_layout <- function(design = NULL) {
    args <- as.list(environment())

    make_interface_function(
        name = "lp_layout",
        args_list = args,
        fn = function(obj) {
            obj@options$layout <- args
            return(obj)
        }
    )
}
