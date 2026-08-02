generate_detected_peaks_for_rects <- function() {
  data.frame(
    mz = c(300.1, 320.4),
    mzmin = c(300.0, 320.3),
    mzmax = c(300.2, 320.5),
    rt = c(2600, 3000),
    rtmin = c(2550, 2950),
    rtmax = c(2650, 3050),
    maxo = c(1e5, 2e5),
    sample_index = c(1, 1),
    sample_id = c("sample_1", "sample_1")
  )
}

rect_layer_data <- function(p) {
  idx <- which(vapply(
    p$layers, function(l) inherits(l$geom, "GeomRect"), logical(1)))
  ggplot2::ggplot_build(p)$data[[idx[1]]]
}

test_that("chrom_peak_rects_layer draws one rectangle per detected peak", {
  # arrange
  peaks <- generate_detected_peaks_for_rects()
  opts <- lcmsPlot:::default_options()
  opts$chrom_peak_rects$show <- TRUE
  base <- ggplot2::ggplot(
    data.frame(rt = c(2500, 3100), mz = c(300, 321)),
    ggplot2::aes(x = rt, y = mz)) + ggplot2::geom_point()

  # act
  p <- lcmsPlot:::chrom_peak_rects_layer(base, peaks, opts)

  # assert
  expect_true(any(vapply(
    p$layers, function(l) inherits(l$geom, "GeomRect"), logical(1))))
  expect_equal(nrow(rect_layer_data(p)), 2L)
})

test_that("chrom_peak_rects_layer keys the rectangles on the peak bounds", {
  # arrange
  peaks <- generate_detected_peaks_for_rects()
  opts <- lcmsPlot:::default_options()
  opts$chrom_peak_rects$show <- TRUE
  base <- ggplot2::ggplot(
    data.frame(rt = c(2500, 3100), mz = c(300, 321)),
    ggplot2::aes(x = rt, y = mz)) + ggplot2::geom_point()

  # act
  rects <- rect_layer_data(
    lcmsPlot:::chrom_peak_rects_layer(base, peaks, opts))

  # assert
  expect_equal(sort(rects$xmin), sort(peaks$rtmin))
  expect_equal(sort(rects$xmax), sort(peaks$rtmax))
})

test_that("chrom_peak_rects_layer follows x_dim when the map is flipped", {
  # arrange
  peaks <- generate_detected_peaks_for_rects()
  opts <- lcmsPlot:::default_options()
  opts$chrom_peak_rects$show <- TRUE
  base <- ggplot2::ggplot(
    data.frame(mz = c(300, 321), rt = c(2500, 3100)),
    ggplot2::aes(x = mz, y = rt)) + ggplot2::geom_point()

  # act
  rects <- rect_layer_data(
    lcmsPlot:::chrom_peak_rects_layer(base, peaks, opts, x_dim = "mz"))

  # assert
  expect_equal(sort(rects$ymin), sort(peaks$rtmin))
  expect_equal(sort(rects$ymax), sort(peaks$rtmax))
})

test_that("chrom_peak_rects_layer gives degenerate peaks a visible height", {
  # arrange: m/z-binned data, where mzmin == mzmax and a plain rect is invisible
  peaks <- generate_detected_peaks_for_rects()
  peaks$mzmin <- peaks$mz
  peaks$mzmax <- peaks$mz
  opts <- lcmsPlot:::default_options()
  opts$chrom_peak_rects$show <- TRUE
  base <- ggplot2::ggplot(
    data.frame(rt = c(2500, 3100), mz = c(300, 321)),
    ggplot2::aes(x = rt, y = mz)) + ggplot2::geom_point()

  # act
  rects <- rect_layer_data(
    lcmsPlot:::chrom_peak_rects_layer(base, peaks, opts))

  # assert
  expect_true(all(rects$ymax > rects$ymin))
})

test_that("chrom_peak_rects_layer clamps the rectangles to the host window", {
  # arrange
  peaks <- generate_detected_peaks_for_rects()
  opts <- lcmsPlot:::default_options()
  opts$chrom_peak_rects$show <- TRUE
  base <- ggplot2::ggplot(
    data.frame(rt = c(2500, 3100), mz = c(300, 321)),
    ggplot2::aes(x = rt, y = mz)) + ggplot2::geom_point()

  # act
  rects <- rect_layer_data(lcmsPlot:::chrom_peak_rects_layer(
    base, peaks, opts, rt_range = c(2600, 3000)))

  # assert
  expect_gte(min(rects$xmin), 2600)
  expect_lte(max(rects$xmax), 3000)
})

test_that("chrom_peak_rects_layer leaves the plot alone when there are no peaks", {
  # arrange
  opts <- lcmsPlot:::default_options()
  opts$chrom_peak_rects$show <- TRUE
  base <- ggplot2::ggplot(
    data.frame(rt = 1, mz = 1), ggplot2::aes(x = rt, y = mz)) +
    ggplot2::geom_point()

  # act / assert
  expect_identical(
    length(lcmsPlot:::chrom_peak_rects_layer(base, NULL, opts)$layers),
    length(base$layers))
  expect_identical(
    length(lcmsPlot:::chrom_peak_rects_layer(
      base, generate_detected_peaks_for_rects()[0, ], opts)$layers),
    length(base$layers))
})

test_that("chrom_peak_rects_layer errors when the m/z bounds are missing", {
  # arrange
  peaks <- generate_detected_peaks_for_rects()
  peaks$mzmin <- NULL
  peaks$mzmax <- NULL
  opts <- lcmsPlot:::default_options()
  opts$chrom_peak_rects$show <- TRUE
  base <- ggplot2::ggplot(
    data.frame(rt = 1, mz = 1), ggplot2::aes(x = rt, y = mz)) +
    ggplot2::geom_point()

  # act / assert
  expect_error(
    lcmsPlot:::chrom_peak_rects_layer(base, peaks, opts),
    "missing the column")
})

test_that("lp_chrom_peak_rects requires an rt / m/z host panel", {
  # arrange
  data_obj <- get_XCMSnExp_object_example(indices = 1:2)

  # act / assert
  expect_error(
    lcmsPlot(data_obj, sample_id_column = "sample_name") +
      lp_chrom_peak_rects(),
    "must be called after")
})

test_that("lp_chrom_peak_rects follows the host panel's samples", {
  # arrange
  data_obj <- get_XCMSnExp_object_example(indices = 1:2, should_group_peaks = TRUE)
  host_samples <- lcmsPlot:::get_metadata(
    data_obj, "sample_name", NULL)$sample_id[1]

  # act: the host restricts the samples, the overlay is left at its default
  obj <- lcmsPlot(data_obj, sample_id_column = "sample_name") +
    lp_intensity_map(
      sample_ids = host_samples,
      mz_range = c(300, 320),
      rt_range = c(2500, 3500)) +
    lp_chrom_peak_rects()

  # assert: inheriting the host's samples keeps the overlay from introducing
  # panels the host never draws
  expect_equal(obj@options$chrom_peak_rects$sample_ids, host_samples)

  # an explicit value still wins
  explicit <- lcmsPlot(data_obj, sample_id_column = "sample_name") +
    lp_intensity_map(
      sample_ids = host_samples,
      mz_range = c(300, 320),
      rt_range = c(2500, 3500)) +
    lp_chrom_peak_rects(sample_ids = "everything")
  expect_equal(explicit@options$chrom_peak_rects$sample_ids, "everything")
})
