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

  # one row per MS1 scan per selected sample, from the scan headers
  scans <- raw_spectra(container, 1L, ms_level = 1L)
  expect_equal(tic$intensity[tic$sample_index == 1], scans$totIonCurrent)
})

test_that("create_total_ion_current sums the peaks with initial = FALSE", {
  # arrange
  container <- lcmsPlotDataContainer(get_test_sample_paths())

  opts <- lcmsPlot:::default_options()
  opts$total_ion_current$sample_ids <- "ko15"
  opts$total_ion_current$initial <- FALSE

  # act
  container <- create_total_ion_current(container, opts)

  # assert: netCDF headers are not the sum of the stored peaks
  tic <- summaries(container)$total_ion_current
  scans <- raw_spectra(container, 1L, ms_level = 1L)
  expect_equal(
    tic$intensity,
    vapply(Spectra::intensity(scans), sum, numeric(1)))
  expect_false(isTRUE(all.equal(tic$intensity, scans$totIonCurrent)))
})

test_that("lp_total_ion_current matches Spectra::tic() on an MsExperiment", {
  # arrange
  data_obj <- get_MsExperiment_object()

  # act
  p <- lcmsPlot(data_obj, sample_id_column = "sample_name") +
    lp_total_ion_current()
  tic <- summaries(p)$total_ion_current

  # assert: the TIC devel plotted for experiment objects
  expected <- Spectra::tic(MsExperiment::spectra(data_obj))
  expect_equal(tic$intensity, expected)
  expect_equal(
    tic$sample_index,
    MsExperiment::spectraSampleIndex(data_obj))
})

test_that("lp_total_ion_current validates initial", {
  # act / assert
  expect_error(lp_total_ion_current(initial = NA), "TRUE or FALSE")
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
