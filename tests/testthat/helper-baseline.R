# Baseline scenarios guarding the data structures behind every plot.
#
# Each scenario builds an lcmsPlotClass object. `capture_plot_datasets()` returns
# the tidy datasets that would be handed to plot_data(), and the reference
# output recorded before the data container refactor (issue #52) lives in
# tests/testthat/fixtures/baseline. tests/baseline/make-baseline.R regenerates
# it.

.baseline_cache <- new.env(parent = emptyenv())

# The xcms pipeline is the slowest fixture, so it is built once per session.
cached_xcms_object <- function(should_group_peaks = TRUE) {
    key <- paste0("xcms_", should_group_peaks)
    if (is.null(.baseline_cache[[key]])) {
        .baseline_cache[[key]] <- get_XCMSnExp_object(
            should_detect_peaks = TRUE,
            should_group_peaks = should_group_peaks)
    }
    .baseline_cache[[key]]
}

baseline_feature_windows <- function() {
    tibble::tibble(
        mzmin = c(334.9, 344.0),
        mzmax = c(335.1, 344.2),
        rtmin = c(2700, 3400),
        rtmax = c(2900, 3500)
    )
}

baseline_feature_ids <- function() {
    head(rownames(xcms::featureDefinitions(cached_xcms_object())), 3)
}

# A LipidSearch result whose lipids sit on two faahKO peaks, so that the
# compound-centric extraction path reads real signal.
baseline_lipid_search_source <- function() {
    paths <- get_test_sample_paths()
    lipids <- data.frame(
        "Rej." = c(0, 0),
        "LipidIon" = c("PC(16:0_18:1)+H", "AcCa(14:0)+H"),
        "LipidGroup" = c("PC(16:0_18:1)", "AcCa(14:0)"),
        "Class" = c("PC", "AcCa"),
        "FattyAcid" = c("(16:0/18:1)", "(14:0)"),
        "CalcMz" = c(335.0, 344.1),
        "IonFormula" = c("C42 H83 O8 N1 P1", "C21 H42 O4 N1"),
        check.names = FALSE
    )
    values <- list(
        "Area" = rbind(c(5000, 3000), c(8000, 6000)),
        "Height" = rbind(c(500, 300), c(800, 600)),
        "Rt" = rbind(c(46.6, 46.7), c(57.5, 57.6)),
        "ObsMz" = rbind(c("335.0", "335.0"), c("344.1", "344.1")),
        "Hwhm(L)" = rbind(c(0.5, 0.5), c(0.4, 0.4)),
        "Hwhm(R)" = rbind(c(0.5, 0.5), c(0.4, 0.4)),
        "Grade" = rbind(c("A", "B"), c("A", "A")),
        "mScore" = rbind(c("8", "7"), c("9", "9")),
        "S/N" = rbind(c("20", "15"), c("30", "30"))
    )
    path <- .write_lipid_search_results(
        file_names = basename(paths), lipids = lipids, values = values)
    LipidSearchSource(path, paths)
}

baseline_scenarios <- list(
    raw_bpc = function() {
        lcmsPlot(get_test_sample_paths()) +
            lp_chromatogram(aggregation_fun = "max")
    },
    raw_tic = function() {
        lcmsPlot(get_test_sample_paths()) +
            lp_chromatogram(aggregation_fun = "sum")
    },
    raw_eic_mass_trace = function() {
        lcmsPlot(get_test_sample_paths()) +
            lp_chromatogram(features = baseline_feature_windows()) +
            lp_mass_trace()
    },
    xcms_feature_ids = function() {
        lcmsPlot(cached_xcms_object(), sample_id_column = "sample_name") +
            lp_chromatogram(
                features = baseline_feature_ids(),
                highlight_peaks = TRUE)
    },
    xcms_rt_both = function() {
        lcmsPlot(cached_xcms_object(), sample_id_column = "sample_name") +
            lp_chromatogram(rt_type = "both")
    },
    xcms_spectra_closest_apex = function() {
        lcmsPlot(cached_xcms_object(), sample_id_column = "sample_name") +
            lp_chromatogram(features = baseline_feature_ids()) +
            lp_spectra(mode = "closest_apex")
    },
    xchromatograms = function() {
        xchr <- xcms::chromatogram(
            cached_xcms_object(should_group_peaks = FALSE),
            mz = as.matrix(baseline_feature_windows()[, c("mzmin", "mzmax")]),
            rt = as.matrix(baseline_feature_windows()[, c("rtmin", "rtmax")]))
        lcmsPlot(xchr) + lp_chromatogram()
    },
    xcms_raw_list = function() {
        lcmsPlot(create_xcms_raw_list(get_test_sample_paths(), profstep = 0)) +
            lp_chromatogram()
    },
    lipid_search = function() {
        lcmsPlot(baseline_lipid_search_source()) +
            lp_lipid_search() +
            lp_chromatogram(highlight_peaks = TRUE)
    },
    total_ion_current = function() {
        lcmsPlot(get_test_sample_paths()) + lp_total_ion_current()
    },
    intensity_map = function() {
        lcmsPlot(get_test_sample_paths()) +
            lp_intensity_map(mz_range = c(330, 350), rt_range = c(2600, 3600))
    },
    rt_diff = function() {
        lcmsPlot(cached_xcms_object(), sample_id_column = "sample_name") +
            lp_rt_diff_plot()
    },
    peak_density = function() {
        lcmsPlot(cached_xcms_object(), sample_id_column = "sample_name") +
            lp_peak_density(features = baseline_feature_windows())
    },
    peak_count_image = function() {
        lcmsPlot(cached_xcms_object(), sample_id_column = "sample_name") +
            lp_peak_count_image()
    }
)

# The tidy datasets handed to plot_data(), keyed by dataset name.
capture_plot_datasets <- function(obj) {
    captured <- NULL
    testthat::local_mocked_bindings(
        plot_data = function(datasets, obj) {
            captured <<- datasets
            NULL
        },
        .package = "lcmsPlot"
    )
    .render_plot(obj, additional_datasets = .additional_datasets(obj))
    captured
}

# Columns whose values must survive the refactor. Keys and helper columns
# (metadata_index, feature_metadata_id, sample_path, ...) are implementation
# details and are left out of the comparison.
.baseline_value_columns <- c(
    "sample_id", "feature_id", "name", "rt_type", "rt", "mz", "intensity",
    "reference", "density", "data_type", "n_peaks", "rt_raw", "rt_adj", "diff"
)

normalise_baseline_dataset <- function(df) {
    cols <- intersect(.baseline_value_columns, colnames(df))
    df <- tibble::as_tibble(as.data.frame(df)[, cols, drop = FALSE])
    df <- df[do.call(order, unname(as.list(df))), , drop = FALSE]
    rownames(df) <- NULL
    df
}

baseline_path <- function(name) {
    testthat::test_path("fixtures", "baseline", paste0(name, ".rds"))
}
