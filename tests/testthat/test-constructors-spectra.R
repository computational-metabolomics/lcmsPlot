test_that("create_spectrum_from_scan_index creates a spectrum from a specified scan index", {
  # arrange
  raw_data <- get_raw_data()
  sample_metadata <- data.frame()
  scan_index <- 321

  # act
  spectrum <- create_spectrum_from_scan_index(
    raw_data,
    sample_meatdata,
    scan_index
  )
  close_raw_data(raw_data)

  # assert
  expect_gt(nrow(spectrum), 0)
  expect_equal(colnames(spectrum), c("mz", "intensity", "rt"))
  expect_equal(unique(spectrum$rt), 3002.163)
})

test_that("create_spectrum_from_scan_index creates a spectrum from a metadata column", {
  # arrange
  raw_data <- get_raw_data()
  sample_metadata <- data.frame(sample_name = "test", scan_index_col = 321)
  scan_index <- "scan_index_col"

  # act
  spectrum <- create_spectrum_from_scan_index(
    raw_data,
    sample_metadata,
    scan_index
  )
  close_raw_data(raw_data)

  # assert
  expect_gt(nrow(spectrum), 0)
  expect_equal(colnames(spectrum), c("mz", "intensity", "rt"))
  expect_equal(unique(spectrum$rt), 3002.163)
})

test_that("create_spectrum_from_closest_scan_to_rt creates a spectrum from the scan closest to a specified RT", {
  # arrange
  raw_data <- get_raw_data()

  # act
  spectrum <- create_spectrum_from_closest_scan_to_rt(
    raw_data,
    rt = 3002,
    ms_level = 1
  )
  close_raw_data(raw_data)

  # assert
  expect_gt(nrow(spectrum), 0)
  expect_equal(colnames(spectrum), c("mz", "intensity", "rt"))
  expect_equal(unique(spectrum$rt), 3002.163)
})
