# `.apply_options()` is the only piece of the Options panel that runs
# outside Shiny's reactive context, so it's also the only piece we can
# unit-test in isolation. The strategy: build a base `lcmsPlotClass` from
# the bundled faahKO CDFs, apply the option-list under test, and inspect
# the `@options` slot to confirm the corresponding `lp_*()` layer ran.

raw_files <- dir(
    system.file("cdf", package = "faahKO"),
    full.names = TRUE,
    recursive  = TRUE)[1:2]

skip_if_no_data <- function() {
    if (length(raw_files) < 2) testthat::skip("faahKO CDFs unavailable")
}

test_that(".apply_options is a no-op when opts is NULL", {
    skip_if_no_data()
    obj <- lcmsPlot(raw_files)
    out <- .apply_options(obj, NULL)
    expect_identical(out, obj)
})

test_that(".apply_options is a no-op when every field is NULL/empty", {
    skip_if_no_data()
    obj <- lcmsPlot(raw_files)
    opts <- list(
        facets          = NULL,
        facet_ncol      = NULL,
        free_x          = FALSE,
        free_y          = FALSE,
        arrange_by      = NULL,
        legend_position = NULL,
        title           = NULL,
        legend_title    = NULL)
    out <- .apply_options(obj, opts)
    expect_identical(out, obj)
})

test_that(".apply_options(facets=...) calls lp_facets", {
    skip_if_no_data()
    obj  <- lcmsPlot(raw_files)
    opts <- list(facets = "sample_id", facet_ncol = 2,
                 free_x = TRUE, free_y = FALSE)
    out  <- .apply_options(obj, opts)
    expect_equal(out@options$facets$facets, "sample_id")
    expect_equal(out@options$facets$ncol,   2)
    expect_true (out@options$facets$free_x)
    expect_false(out@options$facets$free_y)
})

test_that(".apply_options(arrange_by=...) calls lp_arrange", {
    skip_if_no_data()
    obj  <- lcmsPlot(raw_files)
    out  <- .apply_options(obj, list(arrange_by = "sample_id"))
    expect_equal(out@options$arrangement$group_by, "sample_id")
})

test_that(".apply_options(legend_position=...) calls lp_legend", {
    skip_if_no_data()
    obj  <- lcmsPlot(raw_files)
    out  <- .apply_options(obj, list(legend_position = "bottom"))
    expect_equal(out@options$legend$position, "bottom")
})

test_that(".apply_options(title/legend_title=...) calls lp_labels", {
    skip_if_no_data()
    obj  <- lcmsPlot(raw_files)
    out  <- .apply_options(obj, list(title = "Foo", legend_title = "Bar"))
    expect_equal(out@options$labels$title,  "Foo")
    expect_equal(out@options$labels$legend, "Bar")
})

test_that(".empty_to_null returns NULL for '' and the value otherwise", {
    expect_null(.empty_to_null(""))
    expect_null(.empty_to_null(NULL))
    expect_null(.empty_to_null(character(0)))
    expect_equal(.empty_to_null("foo"), "foo")
    expect_equal(.empty_to_null(c("a", "b")), c("a", "b"))
})

test_that(".na_to_null returns NULL for NA and the value otherwise", {
    expect_null(.na_to_null(NA_real_))
    expect_null(.na_to_null(NULL))
    expect_equal(.na_to_null(3),  3)
    expect_equal(.na_to_null(0L), 0L)
})
