test_that("plot_chromatogram plots an extracted ion chromatogram", {
  # arrange
  datasets <- generate_datasets_for_plots()
  supporting_datasets <- list(detected_peaks = TRUE)
  opts <- default_options()
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
  opts <- default_options()
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
  opts <- default_options()
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
