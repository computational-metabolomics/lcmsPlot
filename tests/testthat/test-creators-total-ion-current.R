test_that("create_total_ion_current builds TIC from raw files (character input)", {
  # arrange
  container <- lcmsPlotDataContainer(get_test_sample_paths())

  opts <- lcmsPlot:::default_options()
  opts$total_ion_current$show <- TRUE
  opts$total_ion_current$sample_ids <- c("ko15", "wt15")

  # act
  container <- create_total_ion_current(container, opts)

  # assert
  tic <- summaries(container)$total_ion_current
  expect_gt(nrow(tic), 0)
  expect_equal(colnames(tic), c("intensity", "sample_index"))
  expect_true(is.numeric(tic$intensity))
  expect_setequal(unique(tic$sample_index), sampleMetadata(container)$sample_index)

  # one row per MS1 scan per selected sample, summed from the scan's peaks
  scans <- raw_spectra(container, 1L, ms_level = 1L)
  expect_equal(
    tic$intensity[tic$sample_index == 1],
    vapply(Spectra::intensity(scans), sum, numeric(1)))
})

test_that("create_total_ion_current honours the selected samples", {
  # arrange
  container <- lcmsPlotDataContainer(get_test_sample_paths())

  opts <- lcmsPlot:::default_options()
  opts$total_ion_current$sample_ids <- "wt15"

  # act
  tic <- summaries(create_total_ion_current(container, opts))$total_ion_current

  # assert
  expect_equal(unique(tic$sample_index), 2L)
})

test_that("lp_total_ion_current works end-to-end with raw files", {
  # arrange
  data_obj <- get_test_sample_paths()

  # act
  p <- lcmsPlot(data_obj) +
    lp_total_ion_current(type = "boxplot")

  # assert
  expect_s4_class(p, "lcmsPlotClass")
  expect_gt(nrow(summaries(p)$total_ion_current), 0)
  expect_equal(
    colnames(summaries(p)$total_ion_current),
    c("intensity", "sample_index"))

  # building the ggplot succeeds
  expect_no_error(print(p))
})
