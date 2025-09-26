test_that("create_spectra creates standalone spectra (i.e., without chromatogram)", {
  # arrange
  data_obj <- get_test_sample_paths()
  data_container <- create_data_container_from_obj(data_obj, sample_id_column = "sample_id", metadata = NULL)

  opts <- default_options()
  opts$chromatograms$show <- FALSE
  opts$spectra$show <- TRUE
  opts$spectra$sample_ids <- c('ko15', 'wt15')
  opts$spectra$rt <- 2740
  opts$spectra$mode <- "closest"
  opts$spectra$ms_level <- 1

  # act
  data_container <- create_spectra(data_container, opts)

  # assert
  expect_equal(nrow(data_container@chromatograms), 0)
  expect_gt(nrow(data_container@spectra), 0)
  expect_equal(colnames(data_container@spectra), c("mz", "intensity", "rt", "metadata_index", "additional_metadata_index", "reference"))
  expect_equal(round(unique(data_container@spectra$rt)), 2739)
})
