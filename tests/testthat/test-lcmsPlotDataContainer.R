test_that("lcmsPlotDataContainer standardises the sample metadata", {
  # act
  container <- lcmsPlotDataContainer(get_test_sample_paths())
  samples <- sampleMetadata(container)

  # assert
  expect_s3_class(samples, "tbl_df")
  expect_equal(samples$sample_index, c(1L, 2L))
  expect_equal(samples$sample_id, c("ko15", "wt15"))
  expect_equal(basename(samples$sample_path), c("ko15.CDF", "wt15.CDF"))
})

test_that("lcmsPlotDataContainer only reads the raw spectra when needed", {
  # act
  container <- lcmsPlotDataContainer(get_test_sample_paths())
  loaded <- ensure_spectra(container)

  # assert
  expect_length(MsExperiment::spectra(container@experiment), 0)
  expect_gt(length(MsExperiment::spectra(loaded@experiment)), 0)
  expect_equal(
    unique(raw_spectra(loaded)$sample_index), c(1L, 2L))
  expect_true(all(
    basename(Spectra::dataOrigin(raw_spectra(loaded, 2L))) == "wt15.CDF"))
})

test_that("ensure_spectra reports raw files that do not exist", {
  # arrange
  container <- lcmsPlotDataContainer(c("missing-1.mzML", "missing-2.mzML"))

  # act / assert
  expect_error(ensure_spectra(container), "not found")
})

test_that("ensure_spectra asks for MsBackendRawFileReader to read .raw spectra", {
  skip_if(requireNamespace("MsBackendRawFileReader", quietly = TRUE))
  raw_file <- tempfile(fileext = ".raw")
  file.create(raw_file)
  container <- lcmsPlotDataContainer(raw_file)

  # act / assert
  expect_error(ensure_spectra(container), "MsBackendRawFileReader")
})

test_that("lcmsPlotDataContainer stores XCMSnExp input as an XcmsExperiment", {
  # act
  container <- lcmsPlotDataContainer(
    get_XCMSnExp_object(), sample_id_column = "sample_name")

  # assert
  expect_s4_class(container@source, "XcmsExperiment")
  expect_identical(container@experiment, container@source |>
    (\(x) { MsExperiment::sampleData(x) <- MsExperiment::sampleData(container@experiment); x })())
  expect_equal(sampleMetadata(container)$sample_id, c("ko15", "wt15"))
})

test_that("lcmsPlotDataContainer links XcmsRawList spectra to their samples", {
  # arrange
  raw_list <- create_xcms_raw_list(get_test_sample_paths(), profstep = 0)

  # act
  container <- lcmsPlotDataContainer(raw_list)

  # assert: the in-memory spectra need no loading
  sps <- raw_spectra(container)
  expect_equal(
    as.vector(table(sps$sample_index)),
    vapply(raw_list@data, function(x) length(x@scantime), integer(1)))
  expect_identical(ensure_spectra(container), container)
})

test_that("the container rejects a sample index that is not the row number", {
  # arrange
  container <- lcmsPlotDataContainer(get_test_sample_paths())
  MsExperiment::sampleData(container@experiment)$sample_index <- c(2L, 1L)

  # act / assert
  expect_error(validObject(container), "row number")
})

test_that("the container rejects unknown or malformed summaries", {
  # arrange
  container <- lcmsPlotDataContainer(get_test_sample_paths())

  # act / assert
  expect_error(
    .set_summary(container, "unknown", tibble::tibble(x = 1)),
    "Unknown summaries")
  expect_error(
    .set_summary(container, "total_ion_current", tibble::tibble(x = 1)),
    "intensity")
})

test_that("the accessors work on both the container and the plot object", {
  # arrange
  p <- lcmsPlot(get_test_sample_paths()) +
    lp_chromatogram(features = baseline_feature_windows()) +
    lp_total_ion_current()

  # act / assert
  expect_identical(sampleMetadata(p), sampleMetadata(p@data))
  expect_s4_class(chromatograms(p), "Chromatograms")
  expect_length(chromatograms(p), 4)
  expect_s4_class(spectra(p), "Spectra")
  expect_s3_class(detectedPeaks(p), "data.frame")
  expect_named(summaries(p), "total_ion_current")
})

test_that("mass traces do not depend on the order of the layers", {
  # act
  after <- lcmsPlot(get_test_sample_paths()) +
    lp_chromatogram(features = baseline_feature_windows()) +
    lp_mass_trace()
  before <- lcmsPlot(get_test_sample_paths()) +
    lp_mass_trace() +
    lp_chromatogram(features = baseline_feature_windows())
  chromatograms_only <- lcmsPlot(get_test_sample_paths()) +
    lp_chromatogram(features = baseline_feature_windows())

  # assert
  expect_gt(nrow(summaries(after)$mass_traces), 0)
  expect_equal(summaries(before)$mass_traces, summaries(after)$mass_traces)
  expect_null(summaries(chromatograms_only)$mass_traces)
})

test_that("layers replayed by batch iteration rebuild their data", {
  # arrange
  p <- lcmsPlot(get_test_sample_paths(c(1, 2, 7)), batch_size = 2) +
    lp_total_ion_current() +
    lp_chromatogram()

  # act
  seen <- list()
  iterate_plot_batches(p, function(obj) {
    seen[[length(seen) + 1]] <<- list(
      chromatograms = sort(unique(as_plot_df(chromatograms(obj))$sample_index)),
      tic = sort(unique(summaries(obj)$total_ion_current$sample_index)))
  })

  # assert
  expect_equal(seen[[1]]$chromatograms, c(1L, 2L))
  expect_equal(seen[[2]]$chromatograms, 3L)
  expect_equal(seen[[2]]$tic, c(1L, 2L, 3L))
})
