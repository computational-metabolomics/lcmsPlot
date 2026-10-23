# The values handed to the plotting functions must not change when the data
# structures behind them do. Each scenario is compared against the reference
# recorded by tests/baseline/make-baseline.R.

for (scenario in names(baseline_scenarios)) {
    test_that(paste("plot data matches the baseline:", scenario), {
        expected <- readRDS(baseline_path(scenario))
        actual <- lapply(
            capture_plot_datasets(baseline_scenarios[[scenario]]()),
            normalise_baseline_dataset)

        expect_setequal(names(actual), names(expected))
        for (dataset in names(expected)) {
            expect_equal(
                actual[[dataset]],
                expected[[dataset]],
                tolerance = 1e-6,
                label = paste0(scenario, "$", dataset))
        }
    })
}
