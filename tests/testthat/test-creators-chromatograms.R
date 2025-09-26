test_that("create_chromatograms_from_feature_ids creates correct chromatograms from the grouped peak IDs", {
  # arrange
  data_obj <- get_XCMSnExp_object(should_detect_peaks = TRUE, should_group_peaks = TRUE)
  data_container <- create_data_container_from_obj(data_obj, sample_id_column = "sample_name", metadata = NULL)

  opts <- default_options
  opts$chromatograms$features <- c('M205T2785', 'M207T2713')
  opts$chromatograms$sample_ids <- c('ko15', 'wt15')
  opts$chromatograms$ppm <- 5
  opts$chromatograms$rt_tol <- 10

  # act
  data_container <- create_chromatograms_from_feature_ids(data_container, opts)

  # assert
  expect_gt(nrow(data_container@chromatograms), 0)
  expect_equal(round(min(data_container@chromatograms$rt)), 2705)
  expect_gt(nrow(data_container@mass_traces), 0)
  expect_equal(unique(data_container@additional_metadata$metadata_index), c(1, 2))
  expect_equal(unique(data_container@additional_metadata$feature_id), c('M205T2785', 'M207T2713'))
})

test_that("create_chromatograms_from_features creates correct chromatograms from raw features", {
  # arrange
  data_obj <- get_XCMSnExp_object()
  data_container <- create_data_container_from_obj(data_obj, sample_id_column = "sample_name", metadata = NULL)

  opts <- default_options
  opts$chromatograms$features <- rbind(c(mzmin = 334.9, mzmax = 335.1, rtmin = 2700, rtmax = 2750))
  opts$chromatograms$sample_ids <- c('ko15', 'wt15')

  # act
  data_container <- create_chromatograms_from_features(data_container, opts)

  # assert
  expect_gt(nrow(data_container@chromatograms), 0)
  expect_gte(round(min(data_container@chromatograms$rt)), 2700)
  expect_lte(round(max(data_container@chromatograms$rt)), 2750)
  expect_gt(nrow(data_container@mass_traces), 0)
  expect_equal(unique(data_container@additional_metadata$metadata_index), c(1, 2))
})

test_that("create_full_rt_chromatograms creates BPCs or TICs from an xcms object", {
  # arrange
  data_obj <- get_XCMSnExp_object()
  data_container <- create_data_container_from_obj(data_obj, sample_id_column = "sample_name", metadata = NULL)

  opts <- default_options
  opts$chromatograms$sample_ids <- c('ko15', 'wt15')
  opts$chromatograms$aggregation_fun <- 'max'

  # act
  data_container <- create_full_rt_chromatograms(data_container, opts)

  # assert
  expect_gt(nrow(data_container@chromatograms), 0)
  expect_gte(round(min(data_container@chromatograms$rt)), 2500)
  expect_lte(round(max(data_container@chromatograms$rt)), 4500)
  expect_equal(nrow(data_container@mass_traces), 0)
})
