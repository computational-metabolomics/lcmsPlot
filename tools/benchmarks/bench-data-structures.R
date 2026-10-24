# Time how long each layer takes to compute its plot data.
#
# Run from the package root on two branches and compare:
#   Rscript tools/benchmarks/bench-data-structures.R
#
# Timings cover data creation and the conversion into the tidy datasets handed
# to plot_data(); ggplot construction is excluded.

pkgload::load_all(".", helpers = FALSE, quiet = TRUE)
BiocParallel::register(BiocParallel::SerialParam())

for (helper in list.files(
    "tests/testthat", pattern = "^helper.*\\.R$", full.names = TRUE)) {
    source(helper)
}

paths <- dir(
    system.file("cdf", package = "faahKO"),
    full.names = TRUE,
    recursive = TRUE)[1:3]

set.seed(1)
n_windows <- 100
mz <- runif(n_windows, 210, 590)
rt <- runif(n_windows, 2600, 4400)
windows <- tibble::tibble(
    mzmin = mz - 0.05,
    mzmax = mz + 0.05,
    rtmin = rt - 60,
    rtmax = rt + 60
)

xdata <- suppressMessages(get_XCMSnExp_object(
    should_detect_peaks = TRUE, should_group_peaks = FALSE))

# Windows centred on detected peaks of the first sample, so that every window
# yields an apex spectrum in that sample.
apex_peaks <- xcms::chromPeaks(xdata)
apex_peaks <- head(apex_peaks[apex_peaks[, "sample"] == 1, , drop = FALSE], 20)
apex_windows <- tibble::tibble(
    mzmin = apex_peaks[, "mzmin"] - 0.01,
    mzmax = apex_peaks[, "mzmax"] + 0.01,
    rtmin = apex_peaks[, "rtmin"],
    rtmax = apex_peaks[, "rtmax"],
    sample_id = "ko15"
)

cases <- list(
    eic_100_windows = function() {
        lcmsPlot(paths) + lp_chromatogram(features = windows)
    },
    eic_mass_traces = function() {
        lcmsPlot(paths) +
            lp_chromatogram(features = windows[1:20, ]) +
            lp_mass_trace()
    },
    bpc = function() {
        lcmsPlot(paths) + lp_chromatogram(aggregation_fun = "max")
    },
    tic = function() {
        lcmsPlot(paths) + lp_chromatogram(aggregation_fun = "sum")
    },
    spectra_closest_apex = function() {
        lcmsPlot(xdata, sample_id_column = "sample_name") +
            lp_chromatogram(features = apex_windows, sample_ids = "ko15") +
            lp_spectra(mode = "closest_apex")
    },
    intensity_map = function() {
        lcmsPlot(paths) +
            lp_intensity_map(mz_range = c(200, 600), rt_range = c(2500, 4500))
    },
    total_ion_current = function() {
        lcmsPlot(paths) + lp_total_ion_current()
    }
)

time_case <- function(case, reps = 3) {
    elapsed <- vapply(seq_len(reps), function(i) {
        gc()
        tryCatch(
            system.time(suppressWarnings(suppressMessages(
                capture_plot_datasets(case()))))[["elapsed"]],
            error = function(e) {
                message("  failed: ", conditionMessage(e))
                NA_real_
            })
    }, numeric(1))
    stats::median(elapsed)
}

results <- data.frame(
    case = names(cases),
    median_seconds = vapply(cases, time_case, numeric(1)),
    row.names = NULL
)

branch <- tryCatch(
    system("git rev-parse --abbrev-ref HEAD", intern = TRUE),
    error = function(e) "unknown")
cat("Branch:", branch, "\n")
print(results, row.names = FALSE)
