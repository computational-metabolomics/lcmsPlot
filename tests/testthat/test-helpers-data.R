test_that("get_feature_data returns feature data when specifying mz and rt", {
  # arrange
  opts <- default_options()
  opts$chromatograms$ppm <- 5
  opts$chromatograms$rt_tol <- 10
  feature <- c(mz = 334.908, rt = 2710)
  full_rt_range <- NULL

  # act
  feature_data <- get_feature_data(feature, opts, full_rt_range)

  # assert
  expect_equal(names(feature_data), c("feature_id", "mzr", "rtr"))
  expect_equal(feature_data$mzr, c(334.906325, 334.909675))
  expect_equal(feature_data$rtr, c(2700, 2720))
})

test_that("get_feature_data returns feature data when specifying mz min and max and no rt", {
  # arrange
  opts <- default_options()
  feature <- c(mzmin = 334.906325, mzmax = 334.909675)
  full_rt_range <- c(1, 2000)

  # act
  feature_data <- get_feature_data(feature, opts, full_rt_range)

  # assert
  expect_equal(names(feature_data), c("feature_id", "mzr", "rtr"))
  expect_equal(feature_data$mzr, c(334.906325, 334.909675))
  expect_equal(feature_data$rtr, c(1, 2000))
})
