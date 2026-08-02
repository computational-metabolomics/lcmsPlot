get_xchromatograms_example <- function() {
  data_obj <- get_XCMSnExp_object_example(
    indices = 1:3,
    should_group_peaks = TRUE)

  features <- xcms::featureDefinitions(data_obj)
  features <- features[order(-lengths(features$peakidx)), ][1, ]

  chroms <- xcms::chromatogram(
    data_obj,
    mz = cbind(features$mzmin, features$mzmax),
    rt = cbind(features$rtmin - 60, features$rtmax + 60))

  list(chromatograms = chroms, data_obj = data_obj)
}

test_that("lp_peak_density accepts an XChromatograms object", {
  # arrange
  fixture <- get_xchromatograms_example()

  # act
  obj <- lcmsPlot(fixture$chromatograms) +
    lp_peak_density(bw = 30, min_fraction = 0.4)
  pd <- obj@data@peak_density

  # assert
  expect_gt(nrow(pd), 0)
  expect_true(any(pd$data_type == "density"))
})

test_that("lp_peak_density derives the m/z window from the chromatograms", {
  # arrange: no `features` given, and no lp_chromatogram() to inherit from
  fixture <- get_xchromatograms_example()
  mz_window <- MSnbase::mz(fixture$chromatograms)

  # act
  obj <- lcmsPlot(fixture$chromatograms) + lp_peak_density(bw = 30)

  # assert
  expect_equal(unique(obj@data@peak_density$mzmin), as.numeric(mz_window[, 1]))
  expect_equal(unique(obj@data@peak_density$mzmax), as.numeric(mz_window[, 2]))
})

test_that("lp_peak_density accepts a single XChromatogram", {
  # arrange
  fixture <- get_xchromatograms_example()

  # act
  obj <- lcmsPlot(fixture$chromatograms[1, 1]) +
    lp_peak_density(bw = 30, min_fraction = 0.4)

  # assert
  expect_gt(nrow(obj@data@peak_density), 0)
})

test_that("lp_peak_density(simulate = FALSE) draws the stored features", {
  # arrange
  fixture <- get_xchromatograms_example()
  stored <- as.data.frame(
    xcms::featureDefinitions(fixture$chromatograms))

  # act
  obj <- lcmsPlot(fixture$chromatograms) +
    lp_peak_density(bw = 30, simulate = FALSE)
  rects <- obj@data@peak_density[obj@data@peak_density$data_type == "rect", ]

  # assert
  expect_equal(nrow(rects), nrow(stored))
  expect_equal(as.numeric(rects$rtmin), as.numeric(stored$rtmin))
  expect_equal(as.numeric(rects$rtmax), as.numeric(stored$rtmax))
})

test_that("lp_peak_density still works on an XCMSnExp with explicit features", {
  # arrange
  fixture <- get_xchromatograms_example()
  features <- xcms::featureDefinitions(fixture$data_obj)
  features <- features[order(-lengths(features$peakidx)), ][1, ]

  # act
  obj <- lcmsPlot(fixture$data_obj, sample_id_column = "sample_name") +
    lp_peak_density(
      features = data.frame(
        mzmin = features$mzmin, mzmax = features$mzmax,
        rtmin = features$rtmin, rtmax = features$rtmax),
      bw = 30, min_fraction = 0.4)

  # assert
  expect_gt(nrow(obj@data@peak_density), 0)
})

test_that("lp_peak_density rejects objects without chromatographic peaks", {
  # act / assert
  expect_error(
    lcmsPlot(get_test_sample_paths()) +
      lp_peak_density(features = data.frame(mzmin = 300, mzmax = 320)),
    "no chromatographic peaks")
})
