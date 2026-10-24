# The values handed to the plotting functions must not change when the data
# structures behind them do. Each scenario is compared against the reference
# recorded by tools/baseline/make-baseline.R.

for (scenario in names(baseline_scenarios)) {
    test_that(paste("plot data matches the baseline:", scenario), {
        expected <- readRDS(baseline_path(scenario))
        actual <- lapply(
            capture_plot_datasets(baseline_scenarios[[scenario]]()),
            normalise_baseline_dataset)

        expect_setequal(names(actual), names(expected))
        for (dataset in names(expected)) {
            # New annotation columns may appear; the recorded ones must hold
            # the same values. Columns that were entirely missing carried no
            # values to preserve.
            change <- baseline_changes[[scenario]][[dataset]]
            columns <- setdiff(
                names(Filter(function(x) !all(is.na(x)), expected[[dataset]])),
                change$ignore)
            prepare <- function(df) {
                df <- df[, intersect(columns, colnames(df)), drop = FALSE]
                if (isTRUE(change$distinct)) {
                    df <- distinct(df)
                }
                normalise_baseline_dataset(df)
            }

            expect_equal(
                prepare(actual[[dataset]]),
                prepare(expected[[dataset]]),
                tolerance = 1e-6,
                label = paste0(scenario, "$", dataset))
        }
    })
}
