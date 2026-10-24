# Regenerate the plot-data baselines used by test-baseline-regression.R.
#
# Run from the package root:
#   Rscript tools/baseline/make-baseline.R
#
# The baselines were recorded on devel before the data container refactor
# (issue #52); only regenerate them for deliberate changes in plotted values.

pkgload::load_all(".", helpers = FALSE, quiet = TRUE)
BiocParallel::register(BiocParallel::SerialParam())

for (helper in list.files(
    "tests/testthat", pattern = "^helper.*\\.R$", full.names = TRUE)) {
    source(helper)
}

out_dir <- file.path("tests", "testthat", "fixtures", "baseline")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (name in names(baseline_scenarios)) {
    message("Recording ", name)
    datasets <- local(capture_plot_datasets(baseline_scenarios[[name]]()))
    normalised <- lapply(datasets, normalise_baseline_dataset)
    saveRDS(normalised, file.path(out_dir, paste0(name, ".rds")), version = 3)
}
