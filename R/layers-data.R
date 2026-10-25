#' Define the chromatograms to plot
#'
#' The `lp_chromatogram` function allows the generation of different types of
#' chromatograms.
#'
#' @section Summary chromatograms:
#' In this type of chromatogram, the intensities of the spectra from each scan
#' in an LC–MS dataset are aggregated into a single value per scan.
#' To create such chromatograms do not specify the `features` parameter as
#' that will create the chromatograms for the selected features.
#' In this context, the main parameter is `aggregation_fun`, which can take
#' `max` (base peak chromatogram, BPC), `sum` (total ion current, TIC), or
#' `mean` (averaged ion chromatogram).
#'
#' @section Feature chromatograms:
#' A feature is a combination of retention time (RT) and m/z. Feature
#' chromatograms can be created by specifiying the `features` parameter.
#'
#' @param features Specifies which features to generates the chromatogram for.
#' This can be either:
#' a `matrix` with columns `mz` and `rt` (optional);
#' a `matrix` with columns `mzmin`, `mzmax`,
#' `rtmin` (optional), `rtmax` (optional);
#' a `data.frame` with columns `sample_id`, `mz` and `rt` (optional);
#' a `data.frame` with columns `sample_id`, `mzmin`, `mzmax`,
#' `rtmin` (optional), `rtmax` (optional);
#' a `character` vector of feature identifiers, i.e. the row names of
#' `xcms::featureDefinitions()` - whatever they have been set to, or the
#' `"FT0001"`-style defaults assigned during correspondence. This requires
#' the data to be an `XCMSnExp` or `MsExperiment` object with grouped peaks;
#' identifiers matching no feature are skipped with a warning.
#'
#' For the `matrix` and `data.frame` forms, each feature's identifier (used for
#' labelling and as the `feature_id` faceting column) is taken from a
#' `feature_id` column, or else from the row names (e.g.
#' `rbind(my_feature = c(mz = 335, rt = 2800))`). Features without one get an
#' `xcms::groupnames()`-style `"M<mz>T<rt>"` identifier.
#' @param sample_ids A `character` vector specifying the sample IDs
#' to include in the plot. If `NULL`, the function uses the sample IDs
#' specified in the `lcmsPlot` object.
#' @param ppm A `numeric` value specifying the mass accuracy (in ppm) used
#' when generating chromatograms. Ignored when the `features` parameter
#' specifies both `mzmin` and `mzmax`.
#' @param rt_tol A `numeric` value specifying the RT tolerance used
#' when generating chromatograms. Ignored when the `features` parameter
#' specifies both `rtmin` and `rtmax`.
#' @param line_type A `character` value specifying the line type (from ggplot2).
#' One of: "solid", "dashed", "dotted", "dotdash", "longdash", "twodash".
#' @param highlight_peaks A `logical` value indicating whether to highlight
#' the detected peaks; the input data must be an `XCMSnExp` or `MsExperiment`
#' object.
#' @param highlight_peaks_color A `character` value indicating the color of the
#' highlighted peaks.
#' @param highlight_peaks_mode A `character` value indicating how the peaks
#' should be highlighted. One of `"polygon"` (fills the area under the curve),
#' `"rectangle"` (draws a bounding box from rtmin to rtmax up to the peak apex),
#' or `"point"` (marks the peak apex with a point). Defaults to `"polygon"`.
#' @param highlight_peaks_factor A `character` value indicating the factor from
#' the metadata that determines the color. By default it colors by `sample_id`.
#' @param aggregation_fun A `character` value indicating which aggregation
#' function to use for the spectra intensities; one of `max` (base peak
#' chromatogram), `sum` (total ion current), or `mean` (averaged ion
#' chromatogram). Only applicable to summary chromatograms.
#' @param rt_type A `chracter` value indicating what type of RT to use for the
#' chromatograms. One of `uncorrected` (default), `corrected`, or `both`;
#' the input data must be an `XCMSnExp` or `MsExperiment` object.
#' If `both` is chosen, this will give access to a metadata column called
#' `rt_adjusted` that can be used to differentiate the two RT types (e.g.,
#' through faceting).
#' @param rt_unit A `character` value indicating the unit to use for
#' the RT axis; one of `"minute"` or `"second"`.
#' @param intensity_unit A `character` value indicating the unit to use for
#' the intensity axis; one of `"absolute"` or `"relative"`.
#' @param fill_gaps A `logical` value indicating whether to fill gaps
#' in RT with 0 intensity.
#' @param na.rm A `logical` value. When `TRUE`, data points whose intensity
#' is `NA` are removed before plotting. Defaults to `FALSE`.
#' @param highlight_apices A `logical` value indicating whether to
#' highlight apices with the corresponding RT values in a chromatogram.
#' @param stacked A `numeric` value in `[0, 1]`. When non-zero, each series in a
#' panel is offset up the y axis in m/z order by that fraction of the intensity
#' range, so co-eluting traces stop occluding each other.
#' Defaults to `0` (no stacking).
#' Which series get offset follows `lp_arrange(group_by = )`; without it,
#' samples are stacked.
#' @param transform A `function` applied to the intensities before plotting,
#' such as `log10`, to compress the dynamic range. It is applied to the peak
#' highlight geometry as well, so shaded peaks stay attached to their traces.
#' Defaults to `identity`. Note that `log10` maps the exact zeros injected by
#' `fill_gaps = TRUE` to `-Inf`.
#' @return This function returns another function that takes an `lcmsPlot`
#' object and produces a modified version containing the generated chromatograms
#' in its `data` slot. It is designed to be used with the `+` operator,
#' which serves as a layering mechanism. Each use of `+` incrementally enriches
#' the `lcmsPlot` object by adding new data or visual components.
#' @export
#' @examples
#' raw_files <- dir(
#'    system.file("cdf", package = "faahKO"),
#'    full.names = TRUE,
#'    recursive = TRUE)[1:5]
#'
#' p <- lcmsPlot(raw_files) +
#'   lp_chromatogram(aggregation_fun = "max") +
#'   lp_arrange(group_by = "sample_id")
#'
#' p
lp_chromatogram <- function(
    features = NULL,
    sample_ids = NULL,
    ppm = 10,
    rt_tol = 10,
    line_type = "solid",
    highlight_peaks = FALSE,
    highlight_peaks_color = NULL,
    highlight_peaks_mode = "polygon",
    highlight_peaks_factor = "sample_id",
    aggregation_fun = "max",
    rt_type = "uncorrected",
    rt_unit = "second",
    intensity_unit = "absolute",
    fill_gaps = FALSE,
    na.rm = FALSE,
    highlight_apices = list(column = NULL, top_n = NULL),
    stacked = 0,
    transform = identity
) {
    aggregation_fun <- match.arg(aggregation_fun, c("max", "sum", "mean"))

    if (!is.function(transform)) {
        stop("lp_chromatogram: 'transform' must be a function.")
    }
    if (!is.numeric(stacked) || length(stacked) != 1L) {
        stop("lp_chromatogram: 'stacked' must be a single numeric value.")
    }

    args <- as.list(environment())

    make_interface_function(
        name = "lp_chromatogram",
        args_list = args,
        fn = function(obj) {
            obj@options$chromatograms <- .layer_options(
                args,
                show = TRUE,
                sample_ids = .batch_sample_ids(
                    obj, .default_sample_ids(obj, sample_ids)))

            obj@data <- create_chromatograms(
                obj@data@source,
                obj@data,
                obj@options,
                features)

            if (na.rm) {
                obj@data@chromatograms <- drop_missing_intensities(
                    obj@data@chromatograms)
            }

            # Mass traces are built from the chromatograms' windows, by
            # whichever of the two layers comes second.
            if (isTRUE(obj@options$mass_traces$show)) {
                obj@data <- create_mass_traces(obj@data, obj@options)
            }

            validObject(obj@data)

            return(obj)
        }
    )
}

#' Define the mass trace to plot
#'
#' The `lp_mass_trace` function enables the generation of mass traces,
#' which are graphical representations commonly used in mass spectrometry
#' data analysis. A mass trace plots individual data points defined by their
#' retention time and corresponding mass-to-charge ratio (m/z), making it easier
#' to visualise how specific ions behave over the course of a
#' chromatographic run.
#'
#' @return This function returns another function that takes an `lcmsPlot`
#' object and produces a modified version containing the generated mass traces
#' in its `data` slot. It is designed to be used with the `+` operator,
#' which serves as a layering mechanism. Each use of `+` incrementally enriches
#' the `lcmsPlot` object by adding new data or visual components.
#' @export
#' @examples
#' raw_files <- dir(
#'    system.file("cdf", package = "faahKO"),
#'    full.names = TRUE,
#'    recursive = TRUE)[1:5]
#'
#' ## Create chromatograms of a specific feature
#' p <- lcmsPlot(raw_files) +
#'   lp_chromatogram(features = rbind(c(
#'     mzmin = 334.9,
#'     mzmax = 335.1,
#'     rtmin = 2700,
#'     rtmax = 2900))) +
#'   lp_arrange(group_by = "sample_id")
#'
#' ## Add mass traces
#' p <- p + lp_mass_trace()
#'
#' p
lp_mass_trace <- function() {
    make_interface_function(
        name = "lp_mass_trace",
        args_list = list(),
        fn = function(obj) {
            obj@options$mass_traces$show <- TRUE

            if (length(obj@data@chromatograms) > 0) {
                obj@data <- create_mass_traces(obj@data, obj@options)
            }

            return(obj)
        }
    )
}

#' Define the spectra to plot
#'
#' The `lp_spectra` function enables the generation of spectra, which are
#' graphical representations of ions detected at each mass-to-charge ratio (m/z)
#' with their corresponding absolute or relative intensities.
#'
#' @section Spectra associated with chromatograms:
#' A spectrum is obtained from a scan at a specific retention time (RT).
#' Therefore, when plotting a chromatogram together with its associated spectra,
#' it is common to mark the RT with a vertical line on the chromatogram
#' to indicate where the spectra were acquired. See the example below on
#' how to generate these types of spectra.
#'
#' @section Standalone spectra:
#' Standalone spectra can also be generated, provided no chromatograms
#' are present (i.e., lp_chromatogram has not been used).
#'
#' @param sample_ids A `character` vector specifying the sample IDs
#' to include in the plot. If `NULL`, the function uses the sample IDs
#' specified in the `lcmsPlot` object or the `lp_chromatogram` function.
#' @param mode The method to choose the scan from which to extract the spectra.
#' One of: `closest`, the closest scan to the specified RT - `rt` parameter);
#' `closest_apex`, the closest scan to a detected peak;
#' `across_peak`, selects scans across a detected peak at a certain interval
#' specified in the `interval` parameter.
#' `mode` is not applicable to standalone spectra.
#' @param ms_level The MS level to consider for the scan.
#' @param rt When `mode = "closest"`, the RT to consider.
#' @param scan_index The exact scan index to consider for extracting a spectrum.
#' `scan_index` and `mode` are mutually exclusive.
#' @param interval When `mode = "across_peak."` The RT interval to consider.
#' @param spectral_match_db The database containing reference spectra used
#' for matching and comparison with the input spectra.
#' @param match_target_index The index, ranked by descending match score,
#' identifying which reference spectrum to display in the mirror plot.
#' @param peak_label_size A `numeric` value controlling the font size of
#' m/z labels annotated on spectral peaks. Defaults to `3`.
#' @param intensity_breaks_by A `numeric` value specifying the step size
#' (in percent) between y-axis intensity breaks. Defaults to `20`.
#' @param mz_breaks_n A `numeric` value specifying the approximate number of
#' m/z axis breaks, passed to `scales::pretty_breaks(n = ...)`.
#' Defaults to `6`.
#' @param auto_facet A `logical` value. When `TRUE` (default), a facet is
#' automatically added to separate spectra from different samples or scans.
#' Set to `FALSE` to suppress automatic faceting.
#' @return This function returns another function that takes an `lcmsPlot`
#' object and produces a modified version containing the generated spectra
#' in its `data` slot. It is designed to be used with the `+` operator,
#' which serves as a layering mechanism. Each use of `+` incrementally enriches
#' the `lcmsPlot` object by adding new data or visual components.
#' @export
#' @examples
#' raw_files <- dir(
#'    system.file("cdf", package = "faahKO"),
#'    full.names = TRUE,
#'    recursive = TRUE)[1]
#'
#' p <- lcmsPlot(raw_files) +
#'   lp_chromatogram(features = rbind(c(
#'     mzmin = 334.9,
#'     mzmax = 335.1,
#'     rtmin = 2700,
#'     rtmax = 2900))) +
#'   lp_spectra(mode = "closest", rt = 2785)
#' p
lp_spectra <- function(
    sample_ids = NULL,
    mode = 'closest_apex',
    ms_level = 1,
    rt = NULL,
    scan_index = NULL,
    interval = 3,
    spectral_match_db = NULL,
    match_target_index = NULL,
    peak_label_size = 3,
    intensity_breaks_by = 20,
    mz_breaks_n = 6,
    auto_facet = TRUE
) {
    mode <- match.arg(mode, c("closest_apex", "closest", "across_peak"))

    args <- as.list(environment())

    make_interface_function(
        name = "lp_spectra",
        args_list = args,
        fn = function(obj) {
            # With chromatograms, the spectra follow their samples.
            if (obj@options$chromatograms$show) {
                sample_ids <- obj@options$chromatograms$sample_ids
            }

            obj@options$spectra <- .layer_options(
                args,
                show = TRUE,
                sample_ids = .default_sample_ids(obj, sample_ids))

            obj@data <- create_spectra(obj@data, obj@options)
            return(obj)
        }
    )
}

#' Define the total ion current (TIC)
#'
#' The `lp_total_ion_current` generates summary data for the
#' total ion current (TIC) of the selected samples: one value per MS1 scan.
#' It works both with XCMS objects (`XcmsExperiment`, `MsExperiment`,
#' `XCMSnExp`), `XcmsRawList` objects, and with raw files passed directly to
#' `lcmsPlot()` (a `character` vector of `.mzML`, `.mzXML`, `.CDF`, or `.raw`
#' paths).
#'
#' @param sample_ids A `character` vector specifying the sample IDs
#' to include in the plot. If `NULL`, the function uses the sample IDs
#' specified in the `lcmsPlot` object.
#' @param type A `character` value indicating the type of plot;
#' one of `"boxplot"`, `"violin"`, `"jitter"`.
#' @param initial A `logical` value, as in `Spectra::tic()`. With `TRUE`
#' (the default) each scan's TIC is the one recorded in the file's scan
#' header; with `FALSE` it is the sum of the scan's peak intensities. The two
#' can differ widely: netCDF headers report the instrument's total ion
#' current, which includes signal that is not among the stored peaks. TIC
#' chromatograms (`lp_chromatogram(aggregation_fun = "sum")`) sum the peaks of
#' netCDF files, so on such files they match `initial = FALSE`. Scans whose
#' header reports no TIC are always summed from their peaks.
#' @return This function returns another function that takes an `lcmsPlot`
#' object and produces a modified version containing the generated
#' total ion current (TIC) in its `data` slot. It is designed to be used
#' with the `+` operator, which serves as a layering mechanism.
#' Each use of `+` incrementally enriches the `lcmsPlot` object by
#' adding new data or visual components.
#' @export
#' @examples
#' data_obj <- get_XCMSnExp_object_example()
#'
#' p <- lcmsPlot(data_obj, sample_id_column = "sample_name") +
#'   lp_total_ion_current(type = "violin") +
#'   lp_arrange(group_by = "sample_id")
#' p
lp_total_ion_current <- function(
    sample_ids = NULL,
    type = "boxplot",
    initial = TRUE
) {
    if (!is.logical(initial) || length(initial) != 1L || is.na(initial)) {
        stop("lp_total_ion_current: 'initial' must be TRUE or FALSE.")
    }

    args <- as.list(environment())

    make_interface_function(
        name = "lp_total_ion_current",
        args_list = args,
        fn = function(obj) {
            .require_source(
                obj,
                function(source) {
                    is_xcms_data(source) || is.character(source) ||
                        is(source, "XcmsRawList")
                },
                "total_ion_current: to plot the total ion current the data object should be of class XcmsExperiment, MsExperiment, XCMSnExp, XcmsRawList, or a character vector of raw file paths.")

            obj@options$total_ion_current <- .layer_options(
                args,
                show = TRUE,
                sample_ids = .default_sample_ids(obj, sample_ids))

            obj@data <- create_total_ion_current(obj@data, obj@options)

            return(obj)
        }
    )
}

#' Define a 2D intensity map
#'
#' The `lp_intensity_map` function produces an intensity map
#' in which signal intensity is represented at each m/z / RT coordinate.
#'
#' @param mz_range A `numeric` value indicating the m/z range of the map.
#' @param rt_range A `numeric` value indicating the RT range of the map.
#' @param sample_ids A `character` vector specifying the sample IDs
#' to include in the plot. If `NULL`, the function uses the sample IDs
#' specified in the `lcmsPlot` object.
#' @param x_dim A `character` value indicating which dimension to place on the
#' x-axis. One of `"rt"` (default) or `"mz"`.
#' @param y_dim A `character` value indicating which dimension to place on the
#' y-axis. One of `"mz"` (default) or `"rt"`.
#' @param fill_scale A ggplot2 scale object to use for the fill aesthetic
#' (e.g. `scale_fill_gradient(low = "white", high = "red")`). When `NULL`
#' (default), uses `scale_fill_viridis_c` for intensity and
#' `scale_fill_viridis_d` for density plots. Applies to `geom = "tile"` and
#' `geom = "density"`.
#' @param geom A `character` value selecting how the map is drawn. One of
#' `"tile"` (default) for a binned heatmap, `"point"` for a scatter of the
#' individual centroids, or `"density"` for a
#' 2D kernel density estimate. `"point"` skips the binning entirely, so gaps in
#' the mass traces stay visible rather than being implied away.
#' @param point_size A `numeric` value controlling the point size when
#' `geom = "point"`.
#' @param bin_rt A `numeric` value giving the retention-time bin width in
#' seconds used when `geom = "tile"`. Defaults to `0.1`, which is too fine for
#' slow-scanning data and too coarse for fast-scanning data.
#' @param bin_mz A `numeric` value giving the m/z bin width used when
#' `geom = "tile"`. Defaults to `0.1`; high-resolution data usually wants less.
#' @param colour_scale A ggplot2 scale object to use for the colour aesthetic
#' when `geom = "point"`. When `NULL` (default), uses `scale_colour_viridis_c`.
#' @return This function returns another function that takes an
#' `lcmsPlot` object and produces a modified version containing the generated
#' 2D intensity map in its `data` slot. It is designed to be used with the
#' `+` operator, which serves as a layering mechanism.
#' Each use of `+` incrementally enriches the `lcmsPlot` object by
#' adding new data or visual components.
#' @export
#' @examples
#' raw_files <- dir(
#'   system.file("cdf", package = "faahKO"),
#'   full.names = TRUE,
#'   recursive = TRUE)[1]
#'
#' p <- lcmsPlot(raw_files) +
#'   lp_intensity_map(
#'     mz_range = c(200, 600),
#'     rt_range = c(4200, 4500),
#'     geom = "density")
#' p
lp_intensity_map <- function(
    mz_range,
    rt_range,
    sample_ids = NULL,
    x_dim = "rt",
    y_dim = "mz",
    fill_scale = NULL,
    geom = "tile",
    point_size = 0.5,
    bin_rt = 0.1,
    bin_mz = 0.1,
    colour_scale = NULL
) {
    geom <- match.arg(geom, c("tile", "point", "density"))
    x_dim <- match.arg(x_dim, c("rt", "mz"))
    y_dim <- match.arg(y_dim, c("mz", "rt"))

    args <- as.list(environment())

    make_interface_function(
        name = "lp_intensity_map",
        args_list = args,
        fn = function(obj) {
            obj@options$intensity_maps <- .layer_options(
                args,
                show = TRUE,
                sample_ids = .default_sample_ids(obj, sample_ids))

            obj@data <- create_intensity_map(obj@data, obj@options)

            return(obj)
        }
    )
}

#' Generate the retention time difference plot between raw and adjusted datasets
#'
#' The `lp_rt_diff_plot` function generates the data necessary to plot
#' the difference between the raw and retention time adjusted datasets.
#' Only applicable to `XCMSnExp` and `MsExperiment` objects.
#'
#' @return This function returns another function that takes an
#' `lcmsPlot` object and produces a modified version containing the generated
#' retention time differences, between raw and adjusted,  in its `data` slot.
#' It is designed to be used with the `+` operator, which serves as a layering
#' mechanism. Each use of `+` incrementally enriches the `lcmsPlot` object by
#' adding new data or visual components.
#' @export
#' @examples
#' data_obj <- get_XCMSnExp_object_example(
#'   indices = 1:3,
#'   should_group_peaks = TRUE)
#' p <- lcmsPlot(data_obj, sample_id_column = "sample_name") +
#'   lp_rt_diff_plot()
#' p
lp_rt_diff_plot <- function() {
    make_interface_function(
        name = "lp_rt_diff_plot",
        args_list = list(),
        fn = function(obj) {
            .require_source(
                obj,
                is_xcms_data,
                "lp_rt_diff_plot: to plot the RT differences the data object should be either of class XcmsExperiment or XCMSnExp.")

            if (!has_adjusted_rtime(obj@data)) {
                stop("lp_rt_diff_plot: RT alignment was not performed.")
            }

            obj@options$rt_diff <- list(show = TRUE)
            obj@data <- create_rt_diff(obj@data, obj@options)

            return(obj)
        }
    )
}

#' Plot the peak density for one or more m/z features
#'
#' The `lp_peak_density` function produces a peak density plot.
#' For each supplied feature (m/z bin):
#' - The x axis shows retention time.
#' - The y axis shows sample indices (1 to n), positioned within the density
#'   range.
#' - Detected peaks are plotted as coloured points at each sample's position.
#' - The kernel density estimate of peak apex RTs is drawn as a line.
#' - When `min_fraction` is supplied, the density-descent grouping algorithm
#'   is simulated and feature groups that pass the threshold are highlighted
#'   with semi-transparent rectangles.
#'
#' @param features A `matrix` or `data.frame` with columns `mzmin` and `mzmax`
#' (required) and `rtmin`, `rtmax` (optional). Each row defines one m/z bin.
#' When omitted, the features are taken from `lp_chromatogram` if it has
#' already been called on the same object. For `XChromatograms` and
#' `XChromatogram` objects it is redundant: each extracted ion chromatogram
#' already carries its own m/z and retention-time window, which is used instead.
#' @param bw A `numeric` value specifying the kernel density bandwidth in
#' seconds. Passed to [stats::density()]. Defaults to `30`.
#' @param min_fraction A `numeric` value in `[0, 1]`. When supplied, simulates
#' the `PeakDensityParam` grouping algorithm and draws semi-transparent
#' rectangles for feature groups where at least this fraction of samples (per
#' group) contain a peak. Defaults to `NULL` (no simulation).
#' @param min_samples An `integer` specifying the minimum absolute number of
#' samples per group required to define a feature group. Only used when
#' `min_fraction` is set. Defaults to `1L`.
#' @param sample_groups A vector of length equal to the number of samples
#' assigning each sample to a group (as in `PeakDensityParam`). Defaults to
#' `NULL`, which treats all samples as a single group.
#' @param max_features An `integer` specifying the maximum number of feature
#' group rectangles to draw per m/z bin. Defaults to `50L`.
#' @param rt_unit A `character` value indicating the unit for the RT axis;
#' one of `"second"` (default) or `"minute"`.
#' @param simulate A `logical` value mirroring `xcms::plotChromPeakDensity()`.
#' `TRUE` descends the density curve to derive the feature groups the supplied
#' parameters *would* produce, using `min_fraction` (defaulting to `0.5` when
#' unset). `FALSE` instead draws the feature definitions the object already
#' stores, ignoring `min_fraction` and `min_samples`. Defaults to `NULL`, which
#' simulates whenever `min_fraction` is given.
#' @return This function returns another function that takes an `lcmsPlot`
#' object and produces a modified version containing the generated peak density
#' data in its `data` slot. It is designed to be used with the `+` operator,
#' which serves as a layering mechanism.
#' @export
#' @examples
#' data_obj <- get_XCMSnExp_object_example(
#'   indices = 1:3,
#'   should_group_peaks = TRUE)
#' p <- lcmsPlot(data_obj, sample_id_column = "sample_name") +
#'   lp_peak_density(
#'     features = rbind(c(mzmin = 334.9, mzmax = 335.1,
#'                        rtmin = 2700,  rtmax = 2900)),
#'     bw = 30,
#'     min_fraction = 0.5)
#' p
lp_peak_density <- function(
    features = NULL,
    bw = 30,
    min_fraction = NULL,
    min_samples = 1L,
    sample_groups = NULL,
    max_features = 50L,
    rt_unit = "second",
    simulate = NULL
) {
    args <- as.list(environment())

    make_interface_function(
        name = "lp_peak_density",
        args_list = args,
        fn = function(obj) {
            if (is.null(get_detected_peaks(obj@data@source))) {
                stop(
                    "lp_peak_density: the data object reports no ",
                    "chromatographic peaks. Provide an XcmsExperiment, ",
                    "XCMSnExp, XChromatograms, or XChromatogram object ",
                    "with detected peaks."
                )
            }

            if (is.null(features)) {
                features <- obj@options$chromatograms$features
            }
            if (is.null(features)) {
                # Chromatogram objects already define one m/z window per EIC,
                # which makes `features` redundant for them.
                features <- get_feature_windows(obj@data@source)
            }
            if (is.null(features)) {
                stop(
                    "lp_peak_density: 'features' is missing and no ",
                    "chromatogram features are set. Either provide ",
                    "'features' or call lp_chromatogram() first."
                )
            }

            features <- as.data.frame(features)
            if (!all(c("mzmin", "mzmax") %in% names(features))) {
                stop("lp_peak_density: 'features' must have columns 'mzmin' and 'mzmax'.")
            }
            if (!"rtmin" %in% names(features)) features$rtmin <- NA_real_
            if (!"rtmax" %in% names(features)) features$rtmax <- NA_real_

            obj@options$peak_density <- .layer_options(
                args, show = TRUE, features = features)

            obj@data <- create_peak_density(obj@data, obj@options)

            return(obj)
        }
    )
}

#' Plot chromatographic peak counts per retention-time bin
#'
#' `lp_peak_count_image()` produces a common quality-control view:
#' retention-time bins on the x axis, samples on
#' the y axis, and fill showing how many chromatographic peaks each sample
#' yielded in each bin.
#'
#' Samples are ordered by injection order, and bins
#' containing no peaks are drawn as zero rather than dropped.
#'
#' @param bin_size A `numeric` value giving the retention-time bin width in
#' seconds. Defaults to `30`.
#' @param log A `logical` value. When `TRUE`, counts are shown on a `log2`
#' scale, which is the usual way to read the plot when counts are skewed. Empty
#' bins are left blank rather than being drawn as `-Inf`.
#' @param sample_ids A `character` vector of sample IDs to include.
#' `NULL` uses all samples.
#' @param rt_range An optional length-2 `numeric` limiting the binned
#' retention-time range. When `NULL` (default) the full acquisition range of the
#' data object is used, so a short run shows as empty bins at the edge.
#' @param fill_scale A ggplot2 scale object to use for the fill aesthetic. When
#' `NULL` (default), uses `scale_fill_viridis_c`.
#' @return This function returns another function that takes an `lcmsPlot`
#' object and produces a modified version containing the generated peak counts
#' in its `data` slot. It is designed to be used with the `+` operator,
#' which serves as a layering mechanism.
#' @export
#' @examples
#' data_obj <- get_XCMSnExp_object_example(
#'   indices = 1:3,
#'   should_group_peaks = TRUE)
#'
#' p <- lcmsPlot(data_obj, sample_id_column = "sample_name") +
#'   lp_peak_count_image(bin_size = 30)
#' p
lp_peak_count_image <- function(
    bin_size = 30,
    log = FALSE,
    sample_ids = NULL,
    rt_range = NULL,
    fill_scale = NULL
) {
    if (!is.numeric(bin_size) || length(bin_size) != 1L || bin_size <= 0) {
        stop("lp_peak_count_image: 'bin_size' must be a positive number.")
    }

    args <- as.list(environment())

    make_interface_function(
        name = "lp_peak_count_image",
        args_list = args,
        fn = function(obj) {
            if (is.null(get_detected_peaks(obj@data@source))) {
                stop(
                    "lp_peak_count_image: the data object reports no ",
                    "chromatographic peaks. Provide an XcmsExperiment or ",
                    "XCMSnExp object with detected peaks."
                )
            }

            obj@options$peak_count_image <- .layer_options(args, show = TRUE)

            obj@data <- create_peak_count_image(obj@data, obj@options)

            return(obj)
        }
    )
}

#' Plot chromatographic peak boundaries in the rt / m/z plane
#'
#' `lp_chrom_peak_rects()` draws one rectangle per detected chromatographic
#' peak, spanning `rtmin`-`rtmax` by `mzmin`-`mzmax`. It works in two ways,
#' matching the two xcms methods that draw the same primitive.
#'
#' Called on its own it owns a panel, drawing the rectangles on an otherwise
#' empty retention time / m/z frame as `xcms::plotChromPeaks()` does, with one
#' facet per sample. Called after `lp_mass_trace()` or `lp_intensity_map()` it
#' instead decorates that panel, as `xcms::plot(type = "XIC")` does, and follows
#' the host's axes, window and samples. On an intensity-versus-rt panel the
#' equivalent is
#' `lp_chromatogram(highlight_peaks = TRUE, highlight_peaks_mode = "rectangle")`.
#'
#' The data object has to report chromatographic peaks. On m/z-binned data
#' `mzmin` equals `mzmax`, so each box is given a minimum height and reads as a
#' horizontal segment rather than vanishing.
#'
#' @param sample_ids A `character` vector of sample IDs to include. `NULL`
#' (default) follows the host layer when there is one, so the overlay never
#' introduces panels the host does not draw, and otherwise uses every sample.
#' @param border A `character` colour for the rectangle outline.
#' @param fill A `character` colour for the rectangle interior, or `NA`
#' (default) for unfilled boxes.
#' @param alpha A `numeric` value in `[0, 1]` controlling opacity.
#' @param linewidth A `numeric` value controlling the outline width.
#' @param rt_range An optional length-2 `numeric` restricting the retention-time
#' range, equivalent to `xlim` in `xcms::plotChromPeaks()`. Ignored when the
#' layer decorates an `lp_intensity_map()` panel, which supplies its own window.
#' @param mz_range An optional length-2 `numeric` restricting the m/z range,
#' equivalent to `ylim` in `xcms::plotChromPeaks()`. Ignored when the layer
#' decorates an `lp_intensity_map()` panel.
#' @return A layer function for use with the `+` operator.
#' @export
#' @examples
#' data_obj <- get_XCMSnExp_object_example(
#'   indices = 1:2,
#'   should_group_peaks = TRUE)
#'
#' ## On its own the layer owns the panel, one facet per sample.
#' p <- lcmsPlot(data_obj, sample_id_column = "sample_name") +
#'   lp_chrom_peak_rects(rt_range = c(2500, 3500), mz_range = c(300, 320))
#' p
#'
#' ## After a host layer it decorates that panel instead.
#' p <- lcmsPlot(data_obj, sample_id_column = "sample_name") +
#'   lp_intensity_map(mz_range = c(300, 320), rt_range = c(2500, 3500)) +
#'   lp_chrom_peak_rects()
#' p
lp_chrom_peak_rects <- function(
        sample_ids = NULL,
        border = "#c0392b",
        fill = NA,
        alpha = 0.25,
        linewidth = 0.3,
        rt_range = NULL,
        mz_range = NULL
) {
    args <- as.list(environment())

    make_interface_function(
        name = "lp_chrom_peak_rects",
        args_list = args,
        fn = function(obj) {
            has_map <- isTRUE(obj@options$intensity_maps$show)
            has_trace <- isTRUE(obj@options$mass_traces$show)

            if (has_trace && !has_map &&
                nrow(.summary(obj@data, "mass_traces")) == 0) {
                stop(
                    "lp_chrom_peak_rects: the mass trace panel is empty, so ",
                    "there is nothing to overlay. Mass traces are built from ",
                    "the raw files, which XChromatograms and XChromatogram ",
                    "objects do not carry; use lp_intensity_map() as the host ",
                    "panel instead."
                )
            }

            # lp_intensity_map() does not extract peaks, so unlike the
            # chromatogram layers this overlay may have to populate them.
            if (nrow(obj@data@detected_peaks) == 0) {
                peaks <- .source_peaks(obj@data)

                if (nrow(peaks) == 0) {
                    stop(
                        "lp_chrom_peak_rects: the data object has no detected ",
                        "peaks. Provide an XcmsExperiment or XCMSnExp object ",
                        "with chromatographic peaks."
                    )
                }

                obj@data@detected_peaks <- peaks
                validObject(obj@data)
            }

            # Fail here rather than drawing nothing: an empty host panel is
            # dropped before the renderer runs, so a missing boundary column
            # would otherwise surface as a silently absent overlay.
            peaks <- obj@data@detected_peaks
            bounds <- c("rtmin", "rtmax", "mzmin", "mzmax")
            missing_cols <- bounds[vapply(bounds, function(col) {
                !col %in% names(peaks) || all(is.na(peaks[[col]]))
            }, logical(1))]
            if (length(missing_cols) > 0) {
                stop(
                    "lp_chrom_peak_rects: the detected peaks are missing the ",
                    "column(s) ", paste(missing_cols, collapse = ", "),
                    ", so the peak boundaries cannot be drawn."
                )
            }

            # Follow the host panel's samples by default, or the overlay would
            # contribute peaks for samples the host does not draw and faceting
            # would invent empty panels for them. Standalone there is no host to
            # follow, so every sample is drawn.
            if (is.null(sample_ids)) {
                sample_ids <- if (has_map) {
                    obj@options$intensity_maps$sample_ids
                } else if (has_trace) {
                    obj@options$chromatograms$sample_ids
                } else {
                    NULL
                }
            }

            obj@options$chrom_peak_rects <- .layer_options(
                args, show = TRUE, sample_ids = sample_ids)

            return(obj)
        }
    )
}
