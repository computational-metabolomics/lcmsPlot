spectra_options <- function(...) {
  opts <- lcmsPlot:::default_options()
  opts$chromatograms$show <- FALSE
  opts$spectra$show <- TRUE
  opts$spectra <- utils::modifyList(opts$spectra, list(...))
  opts
}

test_that("create_spectra reads the spectrum of a scan index", {
  # arrange
  container <- lcmsPlotDataContainer(get_test_sample_paths(1))
  opts <- spectra_options(scan_index = 321)

  # act
  sps <- create_spectra(container, opts)@spectra

  # assert
  expect_length(sps, 1)
  expect_gt(length(Spectra::mz(sps)[[1]]), 0)
  expect_equal(round(Spectra::rtime(sps), 3), 3002.163)
})

test_that("create_spectra reads the scan index from a metadata column", {
  # arrange
  container <- lcmsPlotDataContainer(
    get_test_sample_paths(c(1, 7)),
    sample_id_column = "sample_name",
    metadata = data.frame(
      sample_name = c("a", "b"), scan_index_col = c(321, 322)))
  opts <- spectra_options(scan_index = "scan_index_col")

  # act
  sps <- create_spectra(container, opts)@spectra

  # assert
  expect_length(sps, 2)
  expect_equal(sps$sample_index, c(1L, 2L))
  expect_equal(round(Spectra::rtime(sps)[1], 3), 3002.163)
})

test_that("create_spectra reads the scan closest to a retention time", {
  # arrange
  container <- lcmsPlotDataContainer(get_test_sample_paths())
  opts <- spectra_options(
    sample_ids = c("ko15", "wt15"), mode = "closest", rt = 2740)

  # act
  container <- create_spectra(container, opts)
  df <- as_plot_df(container@spectra)

  # assert
  expect_length(container@chromatograms, 0)
  expect_length(container@spectra, 2)
  expect_named(
    df, c("mz", "intensity", "rt", "sample_index", "feature_id", "reference"))
  expect_equal(round(unique(df$rt)), 2739)
  expect_false(any(df$reference))
})

test_that("create_spectra accepts one retention time per sample", {
  # arrange
  container <- lcmsPlotDataContainer(get_test_sample_paths())
  opts <- spectra_options(mode = "closest", rt = c(wt15 = 3000, ko15 = 2740))

  # act
  container <- create_spectra(container, opts)
  sps <- container@spectra

  # assert: each sample's own target, matched to its closest scan
  closest_rt <- function(sample_index, target) {
    rts <- Spectra::rtime(raw_spectra(container, sample_index, ms_level = 1L))
    rts[which.min(abs(rts - target))]
  }
  expect_equal(
    Spectra::rtime(sps),
    c(closest_rt(1L, 2740), closest_rt(2L, 3000)))
})

test_that("create_spectra scales each spectrum to a maximum of 100", {
  # arrange
  container <- lcmsPlotDataContainer(get_test_sample_paths())
  opts <- spectra_options(mode = "closest", rt = 2740)

  # act
  sps <- create_spectra(container, opts)@spectra

  # assert
  expect_equal(
    vapply(Spectra::intensity(sps), max, numeric(1)),
    c(100, 100))
})

test_that("create_spectra labels apex spectra with the feature of their peak", {
  # arrange
  obj <- lcmsPlot(cached_xcms_object(), sample_id_column = "sample_name") +
    lp_chromatogram(features = baseline_feature_ids())
  peaks <- detectedPeaks(obj)

  # act
  obj <- obj + lp_spectra(mode = "closest_apex")
  sps <- spectra(obj)

  # assert: one spectrum per peak, each labelled with its peak's feature
  expect_length(sps, nrow(peaks))
  expect_equal(sort(sps$feature_id), sort(peaks$feature_id))
  expect_setequal(sps$feature_id, baseline_feature_ids())
})

test_that("create_spectra reads spectra across each peak", {
  # arrange
  obj <- lcmsPlot(cached_xcms_object(), sample_id_column = "sample_name") +
    lp_chromatogram(features = baseline_feature_ids()[1])
  peaks <- detectedPeaks(obj)
  n_expected <- sum(vapply(seq_len(nrow(peaks)), function(i) {
    length(seq(peaks$rtmin[i], peaks$rtmax[i], by = 5))
  }, integer(1)))

  # act
  sps <- spectra(obj + lp_spectra(mode = "across_peak", interval = 5))

  # assert
  expect_length(sps, n_expected)
  expect_true(all(sps$feature_id == baseline_feature_ids()[1]))
})

test_that("create_spectra requires a retention time in closest mode", {
  # arrange
  container <- lcmsPlotDataContainer(get_test_sample_paths(1))

  # act / assert
  expect_error(
    create_spectra(container, spectra_options(mode = "closest")),
    "requires 'rt'")
})

test_that("match_library_spectra pairs each query with its best match", {
  # arrange
  peaks <- list(c(100, 200, 300), c(150, 250))
  library <- make_spectra(
    S4Vectors::DataFrame(msLevel = 2L, rtime = NA_real_, name = c("x", "y")),
    mz = peaks,
    intensity = list(c(1, 2, 3), c(5, 1)))
  query <- make_spectra(
    S4Vectors::DataFrame(msLevel = 2L, rtime = c(10, 20)),
    mz = rev(peaks),
    intensity = list(c(5, 1), c(1, 2, 3)))
  query$sample_index <- c(1L, 2L)
  query$feature_id <- c("A", "B")

  # act
  matches <- match_library_spectra(query, library)

  # assert
  expect_equal(matches$name, c("y", "x"))
  expect_equal(Spectra::rtime(matches), c(10, 20))
  expect_equal(matches$sample_index, c(1L, 2L))
  expect_equal(matches$feature_id, c("A", "B"))
  expect_true(all(matches$reference))
})
