generate_peak_count_image <- function() {
  # Two samples x three bins, with one empty bin in the second sample.
  data.frame(
    rt = rep(c(15, 45, 75), times = 2),
    n_peaks = c(4, 7, 2, 3, 0, 5),
    metadata_index = rep(c(1, 2), each = 3),
    feature_metadata_id = NA_real_,
    sample_id = rep(c("wt15", "ko15"), each = 3),
    sample_index = rep(c(1, 2), each = 3)
  )
}

test_that("plot_peak_count_image draws a tile per sample and bin", {
  # arrange
  datasets <- list(peak_count_image = generate_peak_count_image())
  opts <- lcmsPlot:::default_options()

  # act
  p <- plot_peak_count_image(datasets, NULL, opts, single = TRUE)

  # assert
  expect_true(any(vapply(
    p$layers, function(l) inherits(l$geom, "GeomTile"), logical(1))))
  expect_equal(nrow(ggplot2::ggplot_build(p)$data[[1]]), 6L)
})

test_that("plot_peak_count_image orders samples by injection order", {
  # arrange: sample_id sorts the other way round to sample_index
  datasets <- list(peak_count_image = generate_peak_count_image())
  opts <- lcmsPlot:::default_options()

  # act
  p <- plot_peak_count_image(datasets, NULL, opts, single = TRUE)

  # assert
  expect_equal(levels(p$data$sample_label), c("wt15", "ko15"))
})

test_that("plot_peak_count_image blanks empty bins on the log scale", {
  # arrange
  datasets <- list(peak_count_image = generate_peak_count_image())
  opts <- lcmsPlot:::default_options()
  opts$peak_count_image$log <- TRUE

  # act
  p <- plot_peak_count_image(datasets, NULL, opts, single = TRUE)

  # assert
  expect_equal(p$labels$fill, "log2(peaks)")
  expect_false(any(is.infinite(p$data$fill_value)))
  expect_equal(sum(is.na(p$data$fill_value)), 1L)
})

test_that("create_peak_count_image counts peaks per bin and keeps empty bins", {
  # arrange
  data_obj <- get_XCMSnExp_object_example(indices = 1:2, should_group_peaks = TRUE)

  # act
  obj <- lcmsPlot(data_obj, sample_id_column = "sample_name") +
    lp_peak_count_image(bin_size = 30)
  counts <- obj@data@peak_count_image

  # assert
  expect_equal(
    colnames(counts),
    c("rt", "n_peaks", "metadata_index", "feature_metadata_id"))
  expect_equal(sum(counts$n_peaks), nrow(xcms::chromPeaks(data_obj)))
  expect_true(any(counts$n_peaks == 0))
  # every sample is present in every bin
  expect_equal(
    nrow(counts),
    length(unique(counts$rt)) * length(unique(counts$metadata_index)))
})

test_that("create_peak_count_image bins the way xcms::plotChromPeakImage does", {
  # arrange
  data_obj <- get_XCMSnExp_object_example(indices = 1:2, should_group_peaks = TRUE)
  bin_size <- 30

  xlim <- c(floor(min(xcms::rtime(data_obj))),
            ceiling(max(xcms::rtime(data_obj))))
  brks <- seq(xlim[1], xlim[2], by = bin_size)
  if (brks[length(brks)] < xlim[2]) {
    brks <- c(brks, brks[length(brks)] + bin_size)
  }
  pks <- xcms::chromPeaks(data_obj, rt = xlim, msLevel = 1L)
  reference <- vapply(
    split(pks[, "rt"], as.factor(as.integer(pks[, "sample"]))),
    function(z) graphics::hist(z, breaks = brks, plot = FALSE)$counts,
    numeric(length(brks) - 1))

  # act
  counts <- (lcmsPlot(data_obj, sample_id_column = "sample_name") +
    lp_peak_count_image(bin_size = bin_size))@data@peak_count_image

  got <- vapply(
    sort(unique(counts$metadata_index)),
    function(s) {
      sub <- counts[counts$metadata_index == s, ]
      sub$n_peaks[order(sub$rt)]
    },
    numeric(length(brks) - 1))

  # assert
  expect_equal(unname(got), unname(reference))
})

test_that("lp_peak_count_image rejects objects without detected peaks", {
  # arrange
  raw_files <- get_test_sample_paths()

  # act / assert
  expect_error(
    lcmsPlot(raw_files) + lp_peak_count_image(),
    "no chromatographic peaks")
  expect_error(lp_peak_count_image(bin_size = 0), "positive number")
})
