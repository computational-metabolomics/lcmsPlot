spectra_options <- function(...) {
  opts <- lcmsPlot:::default_options()
  opts$chromatograms$show <- FALSE
  opts$spectra$show <- TRUE
  opts$spectra <- utils::modifyList(opts$spectra, list(...))
  opts
}

test_that("create_spectra reads the spectrum of a scan index", {
  # arrange
  container <- lcmsPlotDataContainer(get_test_sample_paths(1))
  opts <- spectra_options(scan_index = 321)

  # act
  sps <- create_spectra(container, opts)@spectra

  # assert
  expect_length(sps, 1)
  expect_gt(length(Spectra::mz(sps)[[1]]), 0)
  expect_equal(round(Spectra::rtime(sps), 3), 3002.163)
})

test_that("create_spectra reads the scan index from a metadata column", {
  # arrange
  container <- lcmsPlotDataContainer(
    get_test_sample_paths(c(1, 7)),
    sample_id_column = "sample_name",
    metadata = data.frame(
      sample_name = c("a", "b"), scan_index_col = c(321, 322)))
  opts <- spectra_options(scan_index = "scan_index_col")

  # act
  sps <- create_spectra(container, opts)@spectra

  # assert
  expect_length(sps, 2)
  expect_equal(sps$sample_index, c(1L, 2L))
  expect_equal(round(Spectra::rtime(sps)[1], 3), 3002.163)
})

test_that("create_spectra reads the scan closest to a retention time", {
  # arrange
  container <- lcmsPlotDataContainer(get_test_sample_paths())
  opts <- spectra_options(
    sample_ids = c("ko15", "wt15"), mode = "closest", rt = 2740)

  # act
  container <- create_spectra(container, opts)
  df <- as_plot_df(container@spectra)

  # assert
  expect_length(container@chromatograms, 0)
  expect_length(container@spectra, 2)
  expect_named(
    df, c("mz", "intensity", "rt", "sample_index", "feature_id", "reference"))
  expect_equal(round(unique(df$rt)), 2739)
  expect_false(any(df$reference))
})

test_that("create_spectra accepts one retention time per sample", {
  # arrange
  container <- lcmsPlotDataContainer(get_test_sample_paths())
  opts <- spectra_options(mode = "closest", rt = c(wt15 = 3000, ko15 = 2740))

  # act
  container <- create_spectra(container, opts)
  sps <- container@spectra

  # assert: each sample's own target, matched to its closest scan
  closest_rt <- function(sample_index, target) {
    rts <- Spectra::rtime(raw_spectra(container, sample_index, ms_level = 1L))
    rts[which.min(abs(rts - target))]
  }
  expect_equal(
    Spectra::rtime(sps),
    c(closest_rt(1L, 2740), closest_rt(2L, 3000)))
})

test_that("create_spectra scales each spectrum to a maximum of 100", {
  # arrange
  container <- lcmsPlotDataContainer(get_test_sample_paths())
  opts <- spectra_options(mode = "closest", rt = 2740)

  # act
  sps <- create_spectra(container, opts)@spectra

  # assert
  expect_equal(
    vapply(Spectra::intensity(sps), max, numeric(1)),
    c(100, 100))
})

test_that("create_spectra labels apex spectra with the feature of their peak", {
  # arrange
  obj <- lcmsPlot(cached_xcms_object(), sample_id_column = "sample_name") +
    lp_chromatogram(features = baseline_feature_ids())
  peaks <- detectedPeaks(obj)

  # act
  obj <- obj + lp_spectra(mode = "closest_apex")
  sps <- spectra(obj)

  # assert: one spectrum per peak, each labelled with its peak's feature
  expect_length(sps, nrow(peaks))
  expect_equal(sort(sps$feature_id), sort(peaks$feature_id))
  expect_setequal(sps$feature_id, baseline_feature_ids())
})

test_that("create_spectra reads spectra across each peak", {
  # arrange
  obj <- lcmsPlot(cached_xcms_object(), sample_id_column = "sample_name") +
    lp_chromatogram(features = baseline_feature_ids()[1])
  peaks <- detectedPeaks(obj)
  n_expected <- sum(vapply(seq_len(nrow(peaks)), function(i) {
    length(seq(peaks$rtmin[i], peaks$rtmax[i], by = 5))
  }, integer(1)))

  # act
  sps <- spectra(obj + lp_spectra(mode = "across_peak", interval = 5))

  # assert
  expect_length(sps, n_expected)
  expect_true(all(sps$feature_id == baseline_feature_ids()[1]))
})

test_that("create_spectra reads spectra at peaks on the adjusted retention times", {
  # arrange: wt15 is aligned to ko15, so its peaks are on adjusted RTs
  obj <- lcmsPlot(cached_xcms_object(), sample_id_column = "sample_name") +
    lp_chromatogram(features = baseline_feature_ids()[1])
  peaks <- detectedPeaks(obj) |> dplyr::filter(sample_index == 2L)
  adjusted_rt <- raw_spectra(obj@data, 2L, ms_level = 1L)$rtime_adjusted
  closest_adjusted <- function(rt) {
    adjusted_rt[which.min(abs(adjusted_rt - rt))]
  }

  # act
  apex <- spectra(obj + lp_spectra(mode = "closest_apex"))
  across <- spectra(obj + lp_spectra(mode = "across_peak", interval = 5))

  # assert
  expect_equal(
    Spectra::rtime(apex[apex$sample_index == 2L]),
    vapply(peaks$rt, closest_adjusted, numeric(1)))
  expect_true(all(
    Spectra::rtime(across[across$sample_index == 2L]) %in% adjusted_rt))
})

test_that("create_spectra requires a retention time in closest mode", {
  # arrange
  container <- lcmsPlotDataContainer(get_test_sample_paths(1))

  # act / assert
  expect_error(
    create_spectra(container, spectra_options(mode = "closest")),
    "requires 'rt'")
})

test_that("match_library_spectra pairs each query with its best match", {
  # arrange
  peaks <- list(c(100, 200, 300), c(150, 250))
  library <- make_spectra(
    S4Vectors::DataFrame(msLevel = 2L, rtime = NA_real_, name = c("x", "y")),
    mz = peaks,
    intensity = list(c(1, 2, 3), c(5, 1)))
  query <- make_spectra(
    S4Vectors::DataFrame(msLevel = 2L, rtime = c(10, 20)),
    mz = rev(peaks),
    intensity = list(c(5, 1), c(1, 2, 3)))
  query$sample_index <- c(1L, 2L)
  query$feature_id <- c("A", "B")

  # act
  matches <- match_library_spectra(query, library)

  # assert
  expect_equal(matches$name, c("y", "x"))
  expect_equal(Spectra::rtime(matches), c(10, 20))
  expect_equal(matches$sample_index, c(1L, 2L))
  expect_equal(matches$feature_id, c("A", "B"))
  expect_true(all(matches$reference))
})

# The DDA standards with one feature per standard, given as extraction
# windows: there are no detected peaks, so the windows are used.
dda_features <- function() {
  data.frame(
    sample_id = c("l-proline-MS2", "l-kynurenine-MS2"),
    feature_id = c("proline", "kynurenine"),
    mz = c(116.0705, 209.0918),
    rt = c(426, 365))
}

test_that("precursor mode keeps the MS2 scan with the most intense precursor", {
  # arrange
  obj <- lcmsPlot(get_dda_sample_paths()) +
    lp_chromatogram(features = dda_features(), rt_tol = 15)

  # act
  sps <- spectra(obj + lp_spectra(mode = "precursor"))

  # assert: one MS2 spectrum per feature, labelled with its feature
  expect_length(sps, 2)
  expect_equal(Spectra::msLevel(sps), c(2L, 2L))
  expect_equal(sps$feature_id, c("proline", "kynurenine"))
  expect_equal(round(Spectra::rtime(sps), 2), c(426.06, 367.49))
  expect_equal(sps$precursor_mz, Spectra::precursorMz(sps))
  expect_equal(round(sps$precursor_mz, 2), c(116.07, 209.09))
})

test_that("precursor mode selects the scans to keep", {
  # arrange
  obj <- lcmsPlot(get_dda_sample_paths()) +
    lp_chromatogram(features = dda_features(), rt_tol = 15)

  # act
  closest <- spectra(obj + lp_spectra(
    mode = "precursor", select = "closest_apex"))
  all <- spectra(obj + lp_spectra(mode = "precursor", select = "all"))
  consensus <- spectra(obj + lp_spectra(
    mode = "precursor", select = "consensus"))

  # assert: the window centre is the apex
  expect_equal(round(Spectra::rtime(closest), 2), c(426.06, 363.72))

  expect_length(all, 6)
  expect_equal(all$feature_id, rep(c("proline", "kynurenine"), each = 3))
  expect_false(is.unsorted(Spectra::rtime(all)[1:3]))

  expect_length(consensus, 2)
  expect_equal(consensus$feature_id, c("proline", "kynurenine"))
  expect_equal(
    Spectra::rtime(consensus),
    as.numeric(tapply(Spectra::rtime(all), all$feature_id, median)[
      c("proline", "kynurenine")]))
  expect_true(all(lengths(consensus) <= lengths(all)[c(1, 4)]))
})

test_that("precursor mode matches the detected peaks of an external source", {
  # arrange: proline is detected in its run, kynurenine and its 13C
  # isotopologue in theirs; every feature is extracted from both runs.
  sample_paths <- get_dda_sample_paths()
  path <- .write_mztab(
    locations = paste0("file:///data/", basename(sample_paths)),
    smf = data.frame(
      SMF_ID = c(1, 2, 3),
      exp_mass_to_charge = c(116.0705, 209.0918, 210.0955),
      charge = c(1, 1, 1),
      retention_time_in_seconds = c(426, 365, 365),
      retention_time_in_seconds_start = c(420, 359, 359),
      retention_time_in_seconds_end = c(432, 371, 371),
      "abundance_assay[1]" = c(24050, NA, NA),
      "abundance_assay[2]" = c(NA, 18700, 2100),
      check.names = FALSE))
  obj <- lcmsPlot(MzTabMSource(path, sample_paths)) +
    lp_mztab(rt_extend = 0) +
    lp_chromatogram()

  # act / assert: one warning lists every feature without an MS2 scan
  expect_warning(
    sps <- spectra(obj + lp_spectra(mode = "precursor")),
    paste0(
      "no MS2 scan matches the precursor of 4 feature\\(s\\), skipped: ",
      "2 \\(l-proline-MS2\\), 3 \\(l-proline-MS2\\), ",
      "3 \\(l-kynurenine-MS2\\), 1 \\(l-kynurenine-MS2\\)$"))
  expect_equal(sps$feature_id, c("1", "2"))
  expect_equal(sps$sample_index, c(1L, 2L))
})

test_that("precursor mode reads MS2 scans of an MsExperiment", {
  # arrange
  experiment <- MsExperiment::readMsExperiment(
    spectraFiles = get_dda_sample_paths("l-kynurenine"),
    sampleData = data.frame(sample_name = "kynurenine"))
  obj <- lcmsPlot(experiment, sample_id_column = "sample_name") +
    lp_chromatogram(
      features = rbind(kynurenine = c(mz = 209.0918, rt = 365)),
      rt_tol = 15)

  # act
  sps <- spectra(obj + lp_spectra(mode = "precursor", select = "all"))

  # assert
  expect_length(sps, 3)
  expect_true(all(sps$feature_id == "kynurenine"))
})

test_that("the most intense precursor falls back to the total ion current", {
  # arrange
  sps <- make_spectra(
    S4Vectors::DataFrame(
      msLevel = 2L,
      rtime = c(1, 2, 3),
      precursorMz = 100,
      precursorIntensity = c(10, NA, 5),
      dataOrigin = "a.mzML"),
    mz = list(c(50, 60), c(50, 60), c(50, 60)),
    intensity = list(c(1, 1), c(1, 1), c(5, 5)))

  # act / assert
  expect_equal(.precursor_intensity(sps), c(2, 2, 10))
  expect_equal(.precursor_intensity(sps[c(1, 3)]), c(10, 5))
})

test_that("precursor spectra are mirrored against a spectral library", {
  # arrange: a library holding the most intense kynurenine spectrum
  obj <- lcmsPlot(get_dda_sample_paths()) +
    lp_chromatogram(features = dda_features(), rt_tol = 15)
  library_spectrum <- Spectra::filterMsLevel(
    Spectra::Spectra(get_dda_sample_paths("l-kynurenine")), 2L)[3]
  library_spectrum$name <- "L-Kynurenine"
  library_path <- tempfile(fileext = ".msp")
  Spectra::export(
    library_spectrum, MsBackendMsp::MsBackendMsp(), file = library_path)

  # act
  df <- as_plot_df(spectra(obj + lp_spectra(
    mode = "precursor", spectral_match_db = library_path)))

  # assert: every reference spectrum shares its query's precursor
  references <- df |> dplyr::filter(reference)
  expect_gt(nrow(references), 0)
  queries <- df |>
    dplyr::filter(!reference) |>
    dplyr::distinct(feature_id, precursor_mz)
  expect_equal(
    unique(references$precursor_mz),
    queries$precursor_mz[match(unique(references$feature_id), queries$feature_id)])
})
