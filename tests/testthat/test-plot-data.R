test_that("as_plot_df lays chromatograms out one row per data point", {
  # arrange
  chromatograms <- as_memory_chromatograms(
    tibble::tibble(
      sample_index = c(1L, 2L),
      feature_id = c("A", "B"),
      mzMin = c(1, 2),
      mzMax = c(1.5, 2.5)),
    list(
      data.frame(rtime = c(1, 2, 3), intensity = c(10, 20, 30)),
      data.frame(rtime = c(4, 5), intensity = c(40, 50))))

  # act
  df <- as_plot_df(chromatograms)

  # assert: the extraction window is left out
  expect_named(df, c("rt", "intensity", "sample_index", "feature_id"))
  expect_equal(df$rt, c(1, 2, 3, 4, 5))
  expect_equal(df$feature_id, c("A", "A", "A", "B", "B"))
  expect_equal(df$sample_index, c(1L, 1L, 1L, 2L, 2L))
})

test_that("as_plot_df lays spectra out one row per peak", {
  # arrange
  sps <- make_spectra(
    S4Vectors::DataFrame(msLevel = 1L, rtime = c(10, 20)),
    mz = list(c(100, 200), 300),
    intensity = list(c(1, 2), 3))
  sps$sample_index <- c(1L, 2L)
  sps$feature_id <- c("A", NA)
  sps$reference <- FALSE

  # act
  df <- as_plot_df(sps)

  # assert
  expect_named(
    df, c("mz", "intensity", "rt", "sample_index", "feature_id", "reference"))
  expect_equal(df$mz, c(100, 200, 300))
  expect_equal(df$rt, c(10, 10, 20))
  expect_equal(df$feature_id, c("A", "A", NA))
})
