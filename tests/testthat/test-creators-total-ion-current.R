test_that("create_total_ion_current builds TIC from raw files (character input)", {
  # arrange
  data_obj <- get_test_sample_paths()
  data_container <- create_data_container_from_obj(
    data_obj, sample_id_column = "sample_id", metadata = NULL)

  opts <- lcmsPlot:::default_options()
  opts$total_ion_current$show <- TRUE
  opts$total_ion_current$sample_ids <- c("ko15", "wt15")

  # act
  data_container <- create_total_ion_current(data_container, opts)

  # assert
  tic <- data_container@total_ion_current
  expect_gt(nrow(tic), 0)
  expect_equal(
    colnames(tic),
    c("intensity", "metadata_index", "feature_metadata_id"))
  expect_true(is.numeric(tic$intensity))

  # metadata_index values map onto the metadata sample_index values
  expect_true(all(
    unique(tic$metadata_index) %in% data_container@metadata$sample_index))

  # one row per MS1 scan per selected sample
  raw <- get_raw_data(1)
  hdr <- ms_header(raw)
  n_ms1 <- nrow(hdr[hdr$msLevel == 1, ])
  close_raw_data(raw)
  expect_equal(sum(tic$metadata_index == 1), n_ms1)
})

test_that("lp_total_ion_current works end-to-end with raw files", {
  # arrange
  data_obj <- get_test_sample_paths()

  # act
  p <- lcmsPlot(data_obj) +
    lp_total_ion_current(type = "boxplot")

  # assert
  expect_s4_class(p, "lcmsPlotClass")
  expect_gt(nrow(p@data@total_ion_current), 0)
  expect_equal(
    colnames(p@data@total_ion_current),
    c("intensity", "metadata_index", "feature_metadata_id"))

  # building the ggplot succeeds
  expect_no_error(print(p))
})
