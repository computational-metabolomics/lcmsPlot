# The Okabe-Ito palette (Okabe & Ito, 2008), distinguishable under the common
# forms of colour blindness.
.okabe_ito <- c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#000000")

.theme_presets <- c("minimal", "classic", "publication")

.palette_presets <- c("okabe_ito", "viridis", "ggplot2")

#' Colours of a discrete palette
#'
#' @param palette A palette preset (`"okabe_ito"`, `"viridis"`, or
#' `"ggplot2"`) or an unnamed `character` vector of colours.
#' @param n The number of colours needed.
#' @return A `character` vector of `n` colours. Okabe-Ito falls back to
#' viridis beyond its 8 colours.
#' @keywords internal
.palette_colours <- function(palette, n) {
    if (n == 0) {
        return(character(0))
    }

    if (length(palette) == 1 && palette %in% .palette_presets) {
        return(switch(
            palette,
            okabe_ito = if (n <= length(.okabe_ito)) {
                .okabe_ito[seq_len(n)]
            } else {
                scales::viridis_pal()(n)
            },
            viridis = scales::viridis_pal()(n),
            ggplot2 = scales::hue_pal()(n)))
    }

    if (length(palette) < n) {
        stop(
            "The palette has ", length(palette), " colours but ", n,
            " groups need colouring. Supply at least ", n, " colours, ",
            "or name them after the groups.", call. = FALSE)
    }
    unname(palette[seq_len(n)])
}

# The values of a column in the order ggplot2 assigns them colours, or NULL
# when the column is continuous and so drawn on a continuous scale.
.discrete_levels <- function(x) {
    if (is.factor(x)) {
        return(levels(droplevels(x)))
    }
    if (!is.character(x) && !is.logical(x)) {
        return(NULL)
    }
    as.character(sort(unique(x)))
}

#' Colours of the series of a plot
#'
#' Resolves the palette once over every dataset, so that a series keeps its
#' colour across all the panels of a patchwork layout, and the highlighted
#' peaks of a series take the colour of its trace.
#'
#' @param datasets A named `list` of the datasets to plot.
#' @param detected_peaks A `data.frame` of the detected peaks, or `NULL`.
#' @param options A `list` of plot options.
#' @return A named `character` vector of colours, keyed by series value, or
#' `NULL` when nothing is coloured by series.
#' @keywords internal
.series_values <- function(datasets, detected_peaks, options) {
    columns <- unique(c(
        options$arrangement$group_by,
        .highlight_column(datasets, detected_peaks, options)))
    if (length(columns) == 0) {
        return(NULL)
    }

    frames <- c(unname(datasets), list(detected_peaks))
    values <- unique(unlist(lapply(columns, function(column) {
        present <- Filter(function(df) column %in% colnames(df), frames)
        if (length(present) == 0) {
            return(NULL)
        }
        .discrete_levels(unlist(
            lapply(present, function(df) df[[column]]),
            use.names = FALSE))
    })))
    if (length(values) == 0) {
        return(NULL)
    }

    palette <- options$theme$palette
    if (!is.null(names(palette))) {
        return(palette)
    }
    stats::setNames(.palette_colours(palette, length(values)), values)
}

# The metadata column that colours the highlighted peaks: the explicit factor,
# else the column that colours the traces.
.highlight_factor <- function(options) {
    options$chromatograms$highlight_peaks_factor %||%
        options$arrangement$group_by %||%
        "sample_id"
}

# The column the highlighted peaks are coloured by, as highlight_peaks() will
# resolve it, or NULL when no peaks are highlighted by colour.
.highlight_column <- function(datasets, detected_peaks, options) {
    chrom_opts <- options$chromatograms
    if (!isTRUE(chrom_opts$highlight_peaks) ||
        !is.null(chrom_opts$highlight_peaks_color)) {
        return(NULL)
    }

    # Polygons are cut out of the chromatograms, the other modes are drawn
    # from the peaks themselves.
    source <- if (identical(chrom_opts$highlight_peaks_mode, "polygon")) {
        datasets$chromatograms
    } else {
        detected_peaks
    }
    column <- .highlight_factor(options)
    if (column %in% colnames(source)) column else "sample_id"
}

# The base font size and line width of the theme, with the preset's defaults
# filled in.
.theme_sizes <- function(options) {
    opts <- options$theme
    base_size <- opts$base_size %||%
        if (identical(opts$preset, "publication")) 8 else 11
    list(
        base_size = base_size,
        line_width = opts$line_width %||% (base_size / 22))
}

# The width of the plotted lines (traces, spectrum sticks).
.line_width <- function(options) {
    if (inherits(options$theme$preset, "theme")) {
        return(options$theme$line_width %||% 0.5)
    }
    .theme_sizes(options)$line_width
}

#' Build the base theme of a plot
#'
#' @param options A `list` of plot options; reads `options$theme`.
#' @return A complete `ggplot2` theme.
#' @keywords internal
.base_theme <- function(options) {
    preset <- options$theme$preset %||% "minimal"
    if (inherits(preset, "theme")) {
        return(preset)
    }

    sizes <- .theme_sizes(options)
    args <- list(
        base_size = sizes$base_size,
        base_line_size = sizes$line_width)

    switch(
        preset,
        minimal = do.call(theme_minimal, args),
        classic = do.call(theme_classic, args),
        publication = do.call(theme_classic, args) +
            theme(
                panel.grid.major = element_line(
                    colour = "grey92", linewidth = rel(0.5)),
                panel.grid.minor = element_blank(),
                panel.background = element_rect(fill = "white", colour = NA),
                plot.background = element_rect(fill = "white", colour = NA),
                strip.background = element_blank(),
                strip.text = element_text(face = "bold"),
                legend.key.size = unit(0.8, "lines")),
        stop("Unknown theme preset: '", preset, "'.", call. = FALSE))
}

#' Colour scales of the series of a plot
#'
#' @param options A `list` of plot options, carrying the series colours
#' resolved by [.series_values()] in `options$theme$series_values`.
#' @param aesthetics The aesthetics coloured by series: `"colour"` and
#' optionally `"fill"`. The fill scale carries no legend of its own.
#' @return A `list` of scales, or `NULL` when nothing is coloured by series.
#' @keywords internal
.series_scales <- function(options, aesthetics = "colour") {
    values <- options$theme$series_values
    if (is.null(values)) {
        return(NULL)
    }

    remove_null_elements(list(
        if ("colour" %in% aesthetics) scale_colour_manual(values = values),
        if ("fill" %in% aesthetics) {
            scale_fill_manual(values = values, guide = "none")
        }))
}

# A discrete scale drawn from the palette, for colourings other than the
# series (e.g. reference against query spectra).
.palette_scale <- function(options, aesthetic = "colour") {
    palette <- options$theme$palette
    if (!is.null(names(palette))) {
        return(scale_colour_manual(values = palette, aesthetics = aesthetic))
    }
    discrete_scale(
        aesthetic,
        palette = function(n) .palette_colours(palette, n))
}
