test_that("create_chromatograms with character features creates correct chromatograms from the grouped peak IDs", {
  # arrange
  data_obj <- get_XCMSnExp_object(should_detect_peaks = TRUE, should_group_peaks = TRUE)
  data_container <- create_data_container_from_obj(data_obj, sample_id_column = "sample_name", metadata = NULL)

  opts <- lcmsPlot:::default_options()
  opts$chromatograms$features <- c('M205T2785', 'M207T2713')
  opts$chromatograms$sample_ids <- c('ko15', 'wt15')
  opts$chromatograms$ppm <- 5
  opts$chromatograms$rt_tol <- 10

  # act
  result <- create_chromatograms(
    data_container@data_obj,
    data_container@metadata,
    opts,
    opts$chromatograms$features
  )
  data_container@chromatograms <- result$chromatograms
  data_container@mass_traces <- result$mass_traces
  data_container@feature_metadata <- result$feature_metadata
  data_container@detected_peaks <- result$detected_peaks

  # assert
  expect_gt(nrow(data_container@chromatograms), 0)
  expect_equal(round(min(data_container@chromatograms$rt)), 2705)
  expect_gt(nrow(data_container@mass_traces), 0)
  expect_equal(unique(data_container@feature_metadata$metadata_index), c(1, 2))
  expect_equal(unique(data_container@feature_metadata$feature_id), c('M205T2785', 'M207T2713'))
})

test_that("create_chromatograms with matrix features creates correct chromatograms from raw features", {
  # arrange
  data_obj <- get_XCMSnExp_object()
  data_container <- create_data_container_from_obj(data_obj, sample_id_column = "sample_name", metadata = NULL)

  opts <- lcmsPlot:::default_options()
  opts$chromatograms$features <- rbind(c(mzmin = 334.9, mzmax = 335.1, rtmin = 2700, rtmax = 2750))
  opts$chromatograms$sample_ids <- c('ko15', 'wt15')

  # act
  result <- create_chromatograms(
    data_container@data_obj,
    data_container@metadata,
    opts,
    opts$chromatograms$features
  )
  data_container@chromatograms <- result$chromatograms
  data_container@mass_traces <- result$mass_traces
  data_container@feature_metadata <- result$feature_metadata
  data_container@detected_peaks <- result$detected_peaks

  # assert
  expect_gt(nrow(data_container@chromatograms), 0)
  expect_gte(round(min(data_container@chromatograms$rt)), 2700)
  expect_lte(round(max(data_container@chromatograms$rt)), 2750)
  expect_gt(nrow(data_container@mass_traces), 0)
  expect_equal(unique(data_container@feature_metadata$metadata_index), c(1, 2))
})

test_that("create_chromatograms with XChromatograms assigns feature_id from mz/rt ranges", {
  # arrange – build an XChromatograms with 2 features × 2 samples
  # XChromatograms (vs MChromatograms) is only returned after peak detection
  xdata <- get_XCMSnExp_object(should_detect_peaks = TRUE)
  mz_ranges <- cbind(mzmin = c(334.9, 206.9), mzmax = c(335.1, 207.1))
  xchrom <- xcms::chromatogram(xdata, mz = mz_ranges)

  data_container <- create_data_container_from_obj(
    xchrom, sample_id_column = NULL, metadata = NULL
  )

  opts <- lcmsPlot:::default_options()
  opts$chromatograms$sample_ids <- data_container@metadata$sample_id

  # act
  result <- create_chromatograms(
    data_container@data_obj,
    data_container@metadata,
    opts,
    NULL
  )

  # assert
  expect_gt(nrow(result$chromatograms), 0)
  # feature_metadata_id must be sequential integers, not NA
  expect_false(any(is.na(result$chromatograms$feature_metadata_id)))
  # feature_metadata must carry a feature_id column
  expect_true("feature_id" %in% colnames(result$feature_metadata))
  # 2 features × 2 samples = 4 rows
  expect_equal(nrow(result$feature_metadata), 4L)
  # feature_id follows M{mz}T{rt} pattern
  expect_match(result$feature_metadata$feature_id[1], "^M\\d+T\\d+$")
  # same feature row → same feature_id across samples
  expect_equal(result$feature_metadata$feature_id[1], result$feature_metadata$feature_id[3])
  expect_equal(result$feature_metadata$feature_id[2], result$feature_metadata$feature_id[4])
  # different rows → different feature_ids
  expect_false(result$feature_metadata$feature_id[1] == result$feature_metadata$feature_id[2])
})

test_that("create_chromatograms with NULL features creates BPCs or TICs from an xcms object", {
  # arrange
  data_obj <- get_XCMSnExp_object()
  data_container <- create_data_container_from_obj(data_obj, sample_id_column = "sample_name", metadata = NULL)

  opts <- lcmsPlot:::default_options()
  opts$chromatograms$sample_ids <- c('ko15', 'wt15')
  opts$chromatograms$aggregation_fun <- 'max'

  # act
  result <- create_chromatograms(
    data_container@data_obj,
    data_container@metadata,
    opts,
    NULL
  )
  data_container@chromatograms <- result$chromatograms
  data_container@mass_traces <- result$mass_traces
  data_container@feature_metadata <- result$feature_metadata
  data_container@detected_peaks <- result$detected_peaks

  # assert
  expect_gt(nrow(data_container@chromatograms), 0)
  expect_gte(round(min(data_container@chromatograms$rt)), 2500)
  expect_lte(round(max(data_container@chromatograms$rt)), 4500)
  expect_equal(nrow(data_container@mass_traces), 0)
})
