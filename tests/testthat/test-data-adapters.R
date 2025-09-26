test_that("get_metadata creates metadata from file paths", {
  # arrange
  paths <- c("test/sample01.mzML", "test/sample02.mzML")
  sample_id_column <- NULL
  metadata <- NULL

  # act
  metadata <- get_metadata(paths, sample_id_column, metadata)

  # assert
  expect_equal(nrow(metadata), 2)
  expect_equal(metadata$sample_index, c(1, 2))
  expect_equal(metadata$sample_id, c("sample01", "sample02"))
})

test_that("get_metadata creates metadata from file paths and existing metadata", {
  # arrange
  paths <- c("test/sample01.mzML", "test/sample02.mzML")
  sample_id_column <- "sample_name"
  metadata <- data.frame(sample_name = c("S1", "S2"))

  # act
  metadata <- get_metadata(paths, sample_id_column, metadata)

  # assert
  expect_equal(nrow(metadata), 2)
  expect_equal(metadata$sample_index, c(1, 2))
  expect_equal(metadata$sample_id, c("S1", "S2"))
  expect_equal(metadata$sample_path, c("test/sample01.mzML", "test/sample02.mzML"))
})

test_that("get_metadata creates metadata from an MsExperiment object", {
  # arrange
  sample_id_column <- "sample_name"
  metadata <- NULL
  obj <- get_MsExperiment_object()

  # act
  metadata <- get_metadata(obj, sample_id_column, metadata)

  # assert
  expect_equal(nrow(metadata), 2)
  expect_equal(metadata$sample_index, c(1, 2))
  expect_equal(metadata$sample_id, c("ko15", "wt15"))
  expect_equal(basename(metadata$sample_path), c("ko15.CDF", "wt15.CDF"))
})

test_that("get_metadata creates metadata from an MsExperiment object and existing metadata", {
  # arrange
  sample_id_column <- "sample_name"
  metadata <- S4Vectors::DataFrame(
    sample_name = c("K1", "W1"),
    test_column = c("S1", "S2")
  )
  obj <- get_MsExperiment_object()

  # act
  metadata <- get_metadata(obj, sample_id_column, metadata)

  # assert
  expect_equal(nrow(metadata), 2)
  expect_equal(metadata$sample_index, c(1, 2))
  expect_equal(metadata$sample_id, c("K1", "W1"))
  expect_equal(basename(metadata$sample_path), c("ko15.CDF", "wt15.CDF"))
  expect_equal(metadata$test_column, c("S1", "S2"))
})

test_that("get_metadata creates metadata from an XCMSnExp object", {
  # arrange
  sample_id_column <- "sample_name"
  metadata <- NULL
  obj <- get_XCMSnExp_object()

  # act
  metadata <- get_metadata(obj, sample_id_column, metadata)

  # assert
  expect_equal(nrow(metadata), 2)
  expect_equal(metadata$sample_index, c(1, 2))
  expect_equal(metadata$sample_id, c("ko15", "wt15"))
  expect_equal(basename(metadata$sample_path), c("ko15.CDF", "wt15.CDF"))
})

test_that("get_detected_peaks returns NULL for input paths", {
  # act
  detected_peaks <- get_detected_peaks(c("sample01.mzML", "sample02.mzML"))

  # assert
  expect_null(detected_peaks)
})

test_that("get_detected_peaks returns the detected peaks for an XCMSnExp object", {
  # arrange
  obj <- get_XCMSnExp_object(should_detect_peaks = TRUE)

  # act
  detected_peaks <- get_detected_peaks(obj)

  # assert
  expect_gt(nrow(detected_peaks), 0)
  expect_contains(colnames(detected_peaks), c("sample_index", "rt", "mz"))
})
