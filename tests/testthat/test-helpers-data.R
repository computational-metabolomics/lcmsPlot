test_that("get_feature_data returns feature data when specifying mz and rt", {
  # arrange
  opts <- lcmsPlot:::default_options()
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
  opts <- lcmsPlot:::default_options()
  feature <- c(mzmin = 334.906325, mzmax = 334.909675)
  full_rt_range <- c(1, 2000)

  # act
  feature_data <- get_feature_data(feature, opts, full_rt_range)

  # assert
  expect_equal(names(feature_data), c("feature_id", "mzr", "rtr"))
  expect_equal(feature_data$mzr, c(334.906325, 334.909675))
  expect_equal(feature_data$rtr, c(1, 2000))
})

test_that("get_feature_data derives an M<mz>T<rt> feature_id when none is given", {
  opts <- lcmsPlot:::default_options()
  feature <- c(mzmin = 334.9, mzmax = 335.1, rtmin = 2700, rtmax = 2900)

  expect_equal(get_feature_data(feature, opts, NULL)$feature_id, "M335T2800")
  expect_equal(get_feature_data(feature, opts, NULL, feature_id = NA)$feature_id, "M335T2800")
})

test_that("get_feature_data uses the supplied feature_id", {
  opts <- lcmsPlot:::default_options()
  feature <- c(mzmin = 334.9, mzmax = 335.1, rtmin = 2700, rtmax = 2900)

  feature_data <- get_feature_data(feature, opts, NULL, feature_id = "mz335")

  expect_equal(feature_data$feature_id, "mz335")
})

test_that("get_feature_names takes matrix row names, ignoring unnamed rows", {
  features <- rbind(
    mz335 = c(mzmin = 334.9, mzmax = 335.1),
    c(mzmin = 278.9, mzmax = 279.1))

  expect_equal(get_feature_names(features), c("mz335", NA))
})

test_that("get_feature_names ignores missing and numeric row names", {
  expect_null(get_feature_names(rbind(c(mz = 335), c(mz = 279))))
  expect_null(get_feature_names(data.frame(mz = c(335, 279))))
  expect_null(get_feature_names(data.frame(mz = c(335, 279))[2, , drop = FALSE]))
})

test_that("get_feature_names prefers a feature_id column over row names", {
  features <- data.frame(mz = c(335, 279), feature_id = c("a", "b"),
                         row.names = c("x", "y"))

  expect_equal(get_feature_names(features), c("a", "b"))
})

test_that("get_features names features from row names", {
  opts <- lcmsPlot:::default_options()
  opts$chromatograms$features <- rbind(
    mz335 = c(mzmin = 334.9, mzmax = 335.1, rtmin = 2700, rtmax = 2900),
    c(mzmin = 278.9, mzmax = 279.1, rtmin = 2740, rtmax = 2840))

  feats <- get_features(opts, tibble::tibble(sample_id = "s1"))

  expect_equal(vapply(feats, `[[`, character(1), "feature_id"), c("mz335", "M279T2790"))
})

test_that("get_features names per-sample features from a feature_id column", {
  opts <- lcmsPlot:::default_options()
  opts$chromatograms$features <- data.frame(
    sample_id = c("s1", "s2", "s2"),
    mz = c(335, 279, 344),
    rt = c(2800, 2790, 2700),
    feature_id = c("a", "b", "c"))

  feats <- get_features(opts, tibble::tibble(sample_id = "s2"))

  expect_equal(vapply(feats, `[[`, character(1), "feature_id"), c("b", "c"))
})
