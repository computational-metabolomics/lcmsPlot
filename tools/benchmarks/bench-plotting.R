# Plotting pipeline benchmarks.
#
# Times the representative workloads on the faahKO CDF files, splitting each
# into data extraction (the `+` calls), plot construction (`.render_plot()`)
# and drawing (`print()` to a null device).
#
# Usage, from the package root:
#   Rscript tools/benchmarks/bench-plotting.R [--save DIR | --compare DIR]
#
# `--save DIR` writes each workload's data slots to DIR as RDS; `--compare DIR`
# checks the current outputs against them with all.equal(), so a refactoring
# can be shown to leave the extracted data unchanged.

suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 2) sub("^--", "", args[1]) else "none"
golden_dir <- if (length(args) >= 2) args[2] else NULL
if (identical(mode, "save")) dir.create(golden_dir, showWarnings = FALSE)

files <- dir(
    system.file("cdf", package = "faahKO"),
    full.names = TRUE,
    recursive = TRUE)

data_slots <- c(
    "chromatograms", "mass_traces", "feature_metadata", "detected_peaks",
    "spectra", "intensity_maps", "total_ion_current", "peak_density",
    "peak_count_image")

results <- list()

run <- function(name, extract) {
    gc()
    t_extract <- system.time(p <- extract())[["elapsed"]]
    t_construct <- system.time(
        p2 <- .render_plot(p, .additional_datasets(p)))[["elapsed"]]
    grDevices::pdf(NULL)
    t_draw <- system.time(print(p2@plot))[["elapsed"]]
    grDevices::dev.off()

    out <- lapply(data_slots, function(s) slot(p@data, s))
    names(out) <- data_slots

    status <- ""
    if (identical(mode, "save")) {
        saveRDS(out, file.path(golden_dir, paste0(name, ".rds")))
    } else if (identical(mode, "compare")) {
        golden <- readRDS(file.path(golden_dir, paste0(name, ".rds")))
        diffs <- vapply(data_slots, function(s) {
            isTRUE(all.equal(
                as.data.frame(golden[[s]]),
                as.data.frame(out[[s]]),
                check.attributes = FALSE))
        }, logical(1))
        status <- if (all(diffs)) {
            "identical"
        } else {
            paste("DIFFERS:", paste(data_slots[!diffs], collapse = ", "))
        }
    }

    results[[name]] <<- data.frame(
        workload = name,
        extract = t_extract,
        construct = t_construct,
        draw = t_draw,
        rows = sum(vapply(out, NROW, integer(1))),
        check = status)
    message(sprintf(
        "%-22s extract %7.3f  construct %6.3f  draw %6.3f  %s",
        name, t_extract, t_construct, t_draw, status))
}

feats5 <- rbind(
    c(mzmin = 334.9, mzmax = 335.1, rtmin = 2700, rtmax = 2900),
    c(mzmin = 300.1, mzmax = 300.3, rtmin = 3300, rtmax = 3500),
    c(mzmin = 508.0, mzmax = 508.2, rtmin = 3600, rtmax = 3800),
    c(mzmin = 344.0, mzmax = 344.2, rtmin = 3900, rtmax = 4100),
    c(mzmin = 462.1, mzmax = 462.3, rtmin = 4100, rtmax = 4300))
feats3 <- rbind(c(mz = 335), c(mz = 300.2), c(mz = 508.1))

run("W1_raw_5feat", function() {
    lcmsPlot(files[1:6]) +
        lp_chromatogram(features = feats5) +
        lp_arrange(group_by = "sample_id")
})
run("W1f_raw_5feat_facets", function() {
    lcmsPlot(files[1:6]) +
        lp_chromatogram(features = feats5) +
        lp_arrange(group_by = "sample_id") +
        lp_facets("feature_id")
})
run("W2_raw_mz_fullrt", function() {
    lcmsPlot(files[1:6]) +
        lp_chromatogram(features = feats3, ppm = 300, fill_gaps = TRUE) +
        lp_arrange(group_by = "sample_id") +
        lp_mass_trace()
})
run("W4_bpc_12", function() {
    lcmsPlot(files) +
        lp_chromatogram() +
        lp_arrange(group_by = "sample_id")
})
run("W6_intmap_tile", function() {
    lcmsPlot(files[1]) +
        lp_intensity_map(mz_range = c(200, 600), rt_range = c(4200, 4500))
})
run("W6b_intmap_point", function() {
    lcmsPlot(files[1]) +
        lp_intensity_map(
            mz_range = c(200, 600), rt_range = c(2500, 4500), geom = "point")
})

xdata <- suppressMessages(get_XCMSnExp_object_example(
    indices = 1:3,
    should_group_peaks = TRUE))
feature_ids <- rownames(xcms::featureDefinitions(xdata))[1:10]
fd <- as.data.frame(xcms::featureDefinitions(xdata))[1:10, ]
feature_matrix <- cbind(
    mzmin = fd$mzmin - 0.01, mzmax = fd$mzmax + 0.01,
    rtmin = fd$rtmin - 20, rtmax = fd$rtmax + 20)

run("W3_xcms_ids_polygon", function() {
    lcmsPlot(xdata, sample_id_column = "sample_name") +
        lp_chromatogram(features = feature_ids, highlight_peaks = TRUE) +
        lp_arrange(group_by = "sample_id")
})
run("W3b_xcms_matrix", function() {
    lcmsPlot(xdata, sample_id_column = "sample_name") +
        lp_chromatogram(features = feature_matrix, highlight_peaks = TRUE) +
        lp_arrange(group_by = "sample_id")
})
run("W5a_spectra_apex", function() {
    lcmsPlot(xdata, sample_id_column = "sample_name") +
        lp_chromatogram(features = feature_ids[1:3]) +
        lp_spectra(mode = "closest_apex")
})
run("W5b_spectra_across", function() {
    lcmsPlot(xdata, sample_id_column = "sample_name") +
        lp_chromatogram(features = feature_ids[1:3]) +
        lp_spectra(mode = "across_peak", interval = 3)
})
run("W7_peak_density", function() {
    lcmsPlot(xdata, sample_id_column = "sample_name") +
        lp_peak_density(features = feature_matrix[1:5, ], min_fraction = 0.5)
})

# Peak highlighting at a size where the per-peak cost shows.
set.seed(1)
ribbon_dataset <- do.call(rbind, lapply(seq_len(20), function(s) {
    tibble::tibble(
        sample_id = paste0("s", s),
        rt = seq(0, 5000, by = 1),
        intensity_plot = stats::runif(5001))
}))
ribbon_peaks <- tibble::tibble(
    sample_id = sample(paste0("s", seq_len(20)), 400, replace = TRUE),
    rtmin = stats::runif(400, 0, 4900))
ribbon_peaks$rtmax <- ribbon_peaks$rtmin + 30

micro <- bench::mark(
    build_ribbon_data = build_ribbon_data(ribbon_dataset, ribbon_peaks),
    iterations = 3,
    check = FALSE)
message("\nMicro-benchmarks:")
print(micro[, c("expression", "min", "median", "mem_alloc")])

message("\nSummary:")
print(do.call(rbind, results), row.names = FALSE)
