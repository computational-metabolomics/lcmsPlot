test_that("plot_chromatogram plots an extracted ion chromatogram", {
  # arrange
  datasets <- generate_datasets_for_plots()
  supporting_datasets <- list(detected_peaks = TRUE)
  opts <- lcmsPlot:::default_options()
  opts$facets$facets <- "sample_id"

  # act
  chrom_plot <- plot_chromatogram(
    datasets,
    supporting_datasets,
    opts,
    single = TRUE
  )

  # assert
  expect_true(any(vapply(chrom_plot$layers, function(layer) inherits(layer$geom, "GeomLine"), logical(1))))
})

test_that("plot_chromatogram plots an extracted ion chromatogram with facets", {
  # arrange
  datasets <- generate_datasets_for_plots()
  supporting_datasets <- list(detected_peaks = TRUE)
  opts <- lcmsPlot:::default_options()
  opts$facets$facets <- "sample_id"

  # act
  chrom_plot <- plot_chromatogram(
    datasets,
    supporting_datasets,
    opts,
    single = TRUE
  )

  # assert
  expect_equal(names(chrom_plot$facet$params$facets), "sample_id")
})

test_that("plot_chromatogram plots an extracted ion chromatogram with a vertical RT line", {
  # arrange
  datasets <- generate_datasets_for_plots()
  supporting_datasets <- list(detected_peaks = TRUE)
  opts <- lcmsPlot:::default_options()
  opts$facets$facets <- "sample_id"
  opts$rt_lines <- list(
    list(
      intercept = 300,
      line_type = "longdash",
      color = "red"
    )
  )

  # act
  chrom_plot <- plot_chromatogram(
    datasets,
    supporting_datasets,
    opts,
    single = TRUE
  )

  # assert
  expect_equal(names(chrom_plot$layers), c("geom_line", "geom_vline"))
  expect_equal(chrom_plot$layers$geom_vline$data$xintercept, 300)
  expect_equal(chrom_plot$layers$geom_vline$aes_params$colour, "red")
  expect_equal(chrom_plot$layers$geom_vline$aes_params$linetype, "longdash")
})

# One row per (sample, feature) so the grouping behaviour is observable.
generate_multi_feature_chromatograms <- function(n_samples = 3, n_features = 3) {
  rt <- seq(0, 20, by = 0.5)
  do.call(rbind, lapply(seq_len(n_samples), function(s) {
    do.call(rbind, lapply(seq_len(n_features), function(f) {
      data.frame(
        rt = rt,
        intensity = dnorm(rt, mean = 4 * f, sd = 0.6) * 100 * f,
        sample_id = paste0("sample_", s),
        sample_group = if (s == 1) "KO" else "WT",
        feature_id = paste0("F", f),
        feature_mz = 100 * f)
    }))
  }))
}

n_series <- function(opts) {
  datasets <- list(chromatograms = generate_multi_feature_chromatograms())
  p <- plot_chromatogram(datasets, list(detected_peaks = TRUE), opts,
                         single = TRUE)
  length(unique(ggplot2::ggplot_build(p)$data[[1]]$group))
}

test_that("plot_chromatogram groups each series separately", {
  # arrange
  opts <- lcmsPlot:::default_options()
  opts$facets$facets <- "sample_id"

  # act / assert: 3 samples, one line each, whatever the colouring
  expect_equal(as.character(opts$arrangement$group_by), character(0))
  expect_equal(n_series(opts), 3L)

  by_sample <- opts
  by_sample$arrangement$group_by <- "sample_id"
  expect_equal(n_series(by_sample), 3L)

  # a coarse grouping factor must not weld samples into one path
  by_group <- opts
  by_group$arrangement$group_by <- "sample_group"
  expect_equal(n_series(by_group), 3L)

  # colouring by feature must split the line, not join the features
  by_feature <- opts
  by_feature$arrangement$group_by <- "feature_id"
  expect_equal(n_series(by_feature), 9L)
})

test_that("plot_chromatogram applies the intensity transform", {
  # arrange
  datasets <- list(chromatograms = generate_multi_feature_chromatograms())
  opts <- lcmsPlot:::default_options()
  opts$facets$facets <- "sample_id"
  logged <- opts
  logged$chromatograms$transform <- log10

  # act
  plain_y <- ggplot2::ggplot_build(plot_chromatogram(
    datasets, list(detected_peaks = TRUE), opts, single = TRUE))$data[[1]]$y
  log_y <- ggplot2::ggplot_build(plot_chromatogram(
    datasets, list(detected_peaks = TRUE), logged, single = TRUE))$data[[1]]$y

  # assert
  expect_equal(max(log_y, na.rm = TRUE), log10(max(plain_y, na.rm = TRUE)))
})

test_that("plot_chromatogram offsets stacked series the way xcms does", {
  # arrange
  datasets <- list(chromatograms = generate_multi_feature_chromatograms())
  opts <- lcmsPlot:::default_options()
  opts$facets$facets <- "sample_id"
  opts$arrangement$group_by <- "feature_id"
  stacked <- opts
  stacked$chromatograms$stacked <- 0.6

  plain_max <- max(ggplot2::ggplot_build(plot_chromatogram(
    datasets, list(detected_peaks = TRUE), opts,
    single = TRUE))$data[[1]]$y, na.rm = TRUE)

  # act
  p <- plot_chromatogram(datasets, list(detected_peaks = TRUE), stacked,
                         single = TRUE)
  offsets <- attr(p, "stacked_offsets")

  # assert: the band spans [0, stacked * ylim], highest m/z on top
  expect_equal(length(offsets), 3L)
  expect_equal(min(offsets), 0)
  expect_equal(max(offsets), 0.6 * plain_max)
  expect_equal(names(offsets)[which.max(offsets)], "F3")
  # a stacked axis has no single meaning, so it is dropped
  expect_null(p$labels$y)
})

test_that("plot_chromatogram leaves output untouched at the defaults", {
  # arrange
  datasets <- list(chromatograms = generate_multi_feature_chromatograms())
  opts <- lcmsPlot:::default_options()
  opts$facets$facets <- "sample_id"
  explicit <- opts
  explicit$chromatograms$stacked <- 0
  explicit$chromatograms$transform <- identity

  # act
  a <- ggplot2::ggplot_build(plot_chromatogram(
    datasets, list(detected_peaks = TRUE), opts, single = TRUE))$data[[1]]
  b <- ggplot2::ggplot_build(plot_chromatogram(
    datasets, list(detected_peaks = TRUE), explicit, single = TRUE))$data[[1]]

  # assert
  expect_equal(a, b)
  expect_equal(names(a), names(b))
})
