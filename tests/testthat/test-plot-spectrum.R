test_that("plot_spectrum plots a spectrum of a scan in a sample", {
  # arrange
  datasets <- generate_datasets_for_plots()
  supporting_datasets <- NULL
  opts <- lcmsPlot:::default_options()
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
  expect_equal(as.character(spectra_plot$layers$geom_text$mapping), "~round(.data$mz, 4)")
  expect_contains(colnames(spectra_plot$layers$geom_text$data), c("mz", "intensity", "rt", "reference", "sample_id"))
})

test_that("plot_spectrum marks the precursor of MS2 spectra", {
  # arrange
  dataset <- tibble::tibble(
    mz = c(70, 116, 70, 116, 50),
    intensity = c(80, 100, 60, 100, 100),
    rt = c(426, 426, 426, 426, 100),
    reference = c(FALSE, FALSE, TRUE, TRUE, FALSE),
    precursor_mz = c(116.0705, 116.0705, 116.0705, 116.0705, NA),
    sample_id = c("a", "a", "a", "a", "b"))
  opts <- lcmsPlot:::default_options()

  # act
  spectra_plot <- plot_spectrum(list(spectra = dataset), NULL, opts)

  # assert: one black diamond on the m/z axis, for the query spectrum only
  expect_contains(names(spectra_plot$layers), "geom_point")
  marker_layer <- spectra_plot$layers$geom_point
  markers <- marker_layer$data
  expect_equal(nrow(markers), 1)
  expect_equal(markers$precursor_mz, 116.0705)
  expect_equal(markers$marker_y, 0)
  expect_equal(marker_layer$aes_params$shape, 23)
  expect_equal(marker_layer$aes_params$colour, "black")
  expect_equal(marker_layer$aes_params$fill, "black")
  expect_setequal(
    unique(spectra_plot$data$sample_id_rt),
    c("a RT: 426 sec. Precursor m/z: 116.0705", "b RT: 100 sec."))
})

test_that("plot_spectrum draws no precursor marker for MS1 spectra", {
  # arrange
  datasets <- generate_datasets_for_plots()
  opts <- lcmsPlot:::default_options()

  # act
  spectra_plot <- plot_spectrum(datasets, NULL, opts)

  # assert
  expect_false("geom_point" %in% names(spectra_plot$layers))
})

test_that("plot_spectrum draws the precursor marker as configured", {
  # arrange
  # a: the precursor left a peak, 0.9 ppm off its m/z; b: it left none
  dataset <- tibble::tibble(
    mz = c(70, 116.0706, 116.2, 70, 150),
    intensity = c(100, 40, 90, 100, 30),
    rt = c(426, 426, 426, 300, 300),
    reference = FALSE,
    precursor_mz = 116.0705,
    sample_id = c("a", "a", "a", "b", "b"))
  opts <- lcmsPlot:::default_options()
  opts$spectra$precursor_marker_shape <- "triangle down filled"
  opts$spectra$precursor_marker_position <- "top"
  opts$spectra$precursor_marker_color <- "red"

  # act
  spectra_plot <- plot_spectrum(list(spectra = dataset), NULL, opts)

  # assert: on top of the precursor's peak, or on the m/z axis without one
  marker_layer <- spectra_plot$layers$geom_point
  markers <- marker_layer$data[order(marker_layer$data$sample_id), ]
  expect_equal(markers$sample_id, c("a", "b"))
  expect_equal(markers$marker_y, c(40, 0))
  expect_equal(marker_layer$aes_params$shape, "triangle down filled")
  expect_equal(marker_layer$aes_params$colour, "red")
  expect_equal(marker_layer$aes_params$fill, "red")
  expect_no_error(ggplot2::ggplot_build(spectra_plot))
})

test_that("lp_spectra rejects an unknown precursor marker position", {
  expect_error(
    lp_spectra(precursor_marker_position = "middle"),
    "'arg' should be one of")
})
