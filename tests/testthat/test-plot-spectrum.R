test_that("plot_spectrum plots a spectrum of a scan in a sample", {
  # arrange
  datasets <- generate_datasets_for_plots()
  supporting_datasets <- NULL
  opts <- default_options()
  opts$facets$facets <- "sample_id"

  # act
  spectra_plot <- plot_spectrum(
    datasets,
    supporting_datasets,
    opts,
    single = TRUE
  )

  # assert
  expect_true(any(vapply(spectra_plot$layers, function(layer) inherits(layer$geom, "GeomSegment"), logical(1))))
  expect_equal(names(spectra_plot$facet$params$facets), "sample_id")
  expect_contains(names(spectra_plot$layers), c("geom_segment", "geom_text"))
  expect_equal(as.character(spectra_plot$layers$geom_text$mapping), "~round(mz, 4)")
  expect_contains(colnames(spectra_plot$layers$geom_text$data), c("mz", "intensity", "rt", "reference", "sample_id"))
})
