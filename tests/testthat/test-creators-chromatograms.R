# Run the chromatogram creator for a container and lay its results out as
# data frames.
run_create_chromatograms <- function(container, opts, features) {
  container <- create_chromatograms(
    container@source, container, opts, features)
  chromatograms <- container@chromatograms
  list(
    container = container,
    chromatograms = as_plot_df(chromatograms),
    chrom_data = if (length(chromatograms) > 0) {
      as.data.frame(Chromatograms::chromData(chromatograms))
    } else {
      data.frame()
    },
    detected_peaks = container@detected_peaks)
}

mass_traces_of <- function(result, opts) {
  summaries(create_mass_traces(result$container, opts))$mass_traces
}

test_that("create_chromatograms with character features creates correct chromatograms from the grouped peak IDs", {
  # arrange
  container <- lcmsPlotDataContainer(
    cached_xcms_object(), sample_id_column = "sample_name")

  opts <- lcmsPlot:::default_options()
  opts$chromatograms$features <- c('FT001', 'FT003')
  opts$chromatograms$sample_ids <- c('ko15', 'wt15')
  opts$chromatograms$ppm <- 5
  opts$chromatograms$rt_tol <- 10

  # act
  result <- run_create_chromatograms(
    container, opts, opts$chromatograms$features)

  # assert
  expect_gt(nrow(result$chromatograms), 0)
  expect_equal(round(min(result$chromatograms$rt)), 2705)
  expect_gt(nrow(mass_traces_of(result, opts)), 0)
  expect_equal(unique(result$chrom_data$sample_index), c(1, 2))
  expect_equal(unique(result$chrom_data$feature_id), c('FT001', 'FT003'))
})

test_that("create_chromatograms with character features extracts one chromatogram per grouped peak", {
  # arrange
  data_obj <- cached_xcms_object()
  container <- lcmsPlotDataContainer(data_obj, sample_id_column = "sample_name")
  feature_ids <- rownames(xcms::featureDefinitions(data_obj))[c(3, 1)]

  opts <- lcmsPlot:::default_options()
  opts$chromatograms$features <- feature_ids
  opts$chromatograms$sample_ids <- sampleMetadata(container)$sample_id
  opts$chromatograms$ppm <- 5
  opts$chromatograms$rt_tol <- 10

  # act
  result <- run_create_chromatograms(
    container, opts, opts$chromatograms$features)

  # assert: the peaks are the features' peakidx rows, features taken in
  # feature definition order
  peak_idx <- xcms::featureDefinitions(data_obj)$peakidx
  names(peak_idx) <- rownames(xcms::featureDefinitions(data_obj))
  feature_ids <- intersect(names(peak_idx), feature_ids)
  expected_rows <- unlist(lapply(feature_ids, function(id) sort(unique(peak_idx[[id]]))))
  all_peaks <- xcms::chromPeaks(data_obj)
  expect_equal(result$detected_peaks$rt, unname(all_peaks[expected_rows, "rt"]))
  expect_equal(result$detected_peaks$sample_index, unname(all_peaks[expected_rows, "sample"]))

  # assert: one chromatogram per peak, drawn from that peak's sample
  expected_features <- rep(
    feature_ids,
    vapply(feature_ids, function(id) length(unique(peak_idx[[id]])), integer(1)))
  expect_equal(nrow(result$chrom_data), length(expected_rows))
  expect_equal(result$chrom_data$sample_index, unname(all_peaks[expected_rows, "sample"]))
  expect_equal(result$chrom_data$feature_id, expected_features)
  expect_equal(result$detected_peaks$feature_id, expected_features)
})

test_that("create_chromatograms with character features honours custom feature definition rownames", {
  # arrange
  data_obj <- cached_xcms_object()
  feat_defs <- xcms::featureDefinitions(data_obj)
  rownames(feat_defs) <- paste0("MyFeature_", seq_len(nrow(feat_defs)))
  xcms::featureDefinitions(data_obj) <- feat_defs

  container <- lcmsPlotDataContainer(data_obj, sample_id_column = "sample_name")

  opts <- lcmsPlot:::default_options()
  opts$chromatograms$features <- c('MyFeature_1', 'MyFeature_3')
  opts$chromatograms$sample_ids <- c('ko15', 'wt15')
  opts$chromatograms$ppm <- 5
  opts$chromatograms$rt_tol <- 10

  # act
  result <- run_create_chromatograms(
    container, opts, opts$chromatograms$features)

  # assert
  expect_gt(nrow(result$chromatograms), 0)
  expect_equal(unique(result$chrom_data$feature_id), c('MyFeature_1', 'MyFeature_3'))
})

test_that("create_chromatograms warns about feature IDs that match no feature", {
  # arrange
  container <- lcmsPlotDataContainer(
    cached_xcms_object(), sample_id_column = "sample_name")

  opts <- lcmsPlot:::default_options()
  opts$chromatograms$features <- c('FT001', 'M205T2785')
  opts$chromatograms$sample_ids <- c('ko15', 'wt15')
  opts$chromatograms$ppm <- 5
  opts$chromatograms$rt_tol <- 10

  # act / assert
  expect_warning(
    result <- run_create_chromatograms(
      container, opts, opts$chromatograms$features),
    "M205T2785"
  )

  # the matching feature is still plotted
  expect_equal(unique(result$chrom_data$feature_id), 'FT001')
})

test_that("create_chromatograms with matrix features creates correct chromatograms from raw features", {
  # arrange
  container <- lcmsPlotDataContainer(
    get_XCMSnExp_object(), sample_id_column = "sample_name")

  opts <- lcmsPlot:::default_options()
  opts$chromatograms$features <- rbind(c(mzmin = 334.9, mzmax = 335.1, rtmin = 2700, rtmax = 2750))
  opts$chromatograms$sample_ids <- c('ko15', 'wt15')

  # act
  result <- run_create_chromatograms(
    container, opts, opts$chromatograms$features)

  # assert
  expect_gt(nrow(result$chromatograms), 0)
  expect_gte(round(min(result$chromatograms$rt)), 2700)
  expect_lte(round(max(result$chromatograms$rt)), 2750)
  expect_gt(nrow(mass_traces_of(result, opts)), 0)
  expect_equal(unique(result$chrom_data$sample_index), c(1, 2))
})

test_that("create_chromatograms with matrix features keeps the peaks inside the windows", {
  # arrange
  data_obj <- cached_xcms_object(should_group_peaks = FALSE)
  container <- lcmsPlotDataContainer(data_obj, sample_id_column = "sample_name")
  window <- c(mzmin = 334.9, mzmax = 335.1, rtmin = 2700, rtmax = 2900)

  opts <- lcmsPlot:::default_options()
  opts$chromatograms$features <- rbind(unname(window))
  colnames(opts$chromatograms$features) <- names(window)
  opts$chromatograms$sample_ids <- c('ko15', 'wt15')

  # act
  result <- run_create_chromatograms(
    container, opts, opts$chromatograms$features)

  # assert
  all_peaks <- as.data.frame(xcms::chromPeaks(data_obj))
  inside <- all_peaks[
    all_peaks$mz >= window[["mzmin"]] & all_peaks$mz <= window[["mzmax"]] &
      all_peaks$rt >= window[["rtmin"]] & all_peaks$rt <= window[["rtmax"]], ]
  expect_gt(nrow(inside), 0)
  expect_setequal(result$detected_peaks$rt, inside$rt)
  expect_true(all(result$detected_peaks$feature_id == "M335T2800"))
})

test_that("create_chromatograms with XChromatograms assigns feature_id from mz/rt ranges", {
  # arrange – build an XChromatograms with 2 features × 2 samples
  # XChromatograms (vs MChromatograms) is only returned after peak detection
  xdata <- cached_xcms_object(should_group_peaks = FALSE)
  mz_ranges <- cbind(mzmin = c(334.9, 206.9), mzmax = c(335.1, 207.1))
  xchrom <- xcms::chromatogram(xdata, mz = mz_ranges)

  container <- lcmsPlotDataContainer(xchrom, sample_id_column = NULL)

  opts <- lcmsPlot:::default_options()
  opts$chromatograms$sample_ids <- sampleMetadata(container)$sample_id

  # act
  result <- run_create_chromatograms(container, opts, NULL)
  feature_id <- result$chrom_data$feature_id

  # assert
  expect_gt(nrow(result$chromatograms), 0)
  # 2 features × 2 samples = 4 chromatograms
  expect_equal(nrow(result$chrom_data), 4L)
  # feature_id follows M{mz}T{rt} pattern
  expect_match(feature_id[1], "^M\\d+T\\d+$")
  # same feature row → same feature_id across samples
  expect_equal(feature_id[1], feature_id[3])
  expect_equal(feature_id[2], feature_id[4])
  # different rows → different feature_ids
  expect_false(feature_id[1] == feature_id[2])
  expect_equal(result$chrom_data$sample_index, c(1L, 1L, 2L, 2L))
})

test_that("create_chromatograms with XChromatograms uses user-set row names", {
  xdata <- cached_xcms_object(should_group_peaks = FALSE)
  mz_ranges <- cbind(mzmin = c(334.9, 206.9), mzmax = c(335.1, 207.1))
  xchrom <- xcms::chromatogram(xdata, mz = mz_ranges)
  rownames(xchrom) <- c("mz335", "mz207")

  container <- lcmsPlotDataContainer(xchrom, sample_id_column = NULL)

  opts <- lcmsPlot:::default_options()
  opts$chromatograms$sample_ids <- sampleMetadata(container)$sample_id

  result <- run_create_chromatograms(container, opts, NULL)

  expect_equal(result$chrom_data$feature_id, rep(c("mz335", "mz207"), times = 2))
})

test_that("create_chromatograms with featureChromatograms() output uses the feature definition row names", {
  xdata <- cached_xcms_object()
  feature_ids <- rownames(xcms::featureDefinitions(xdata))[c(3, 5)]
  xchrom <- xcms::featureChromatograms(xdata, features = feature_ids)

  container <- lcmsPlotDataContainer(xchrom, sample_id_column = NULL)

  opts <- lcmsPlot:::default_options()
  opts$chromatograms$sample_ids <- sampleMetadata(container)$sample_id

  result <- run_create_chromatograms(container, opts, NULL)

  expect_equal(result$chrom_data$feature_id, rep(feature_ids, times = 2))
})

test_that("create_chromatograms with NULL features creates BPCs or TICs from an xcms object", {
  # arrange
  container <- lcmsPlotDataContainer(
    get_XCMSnExp_object(), sample_id_column = "sample_name")

  opts <- lcmsPlot:::default_options()
  opts$chromatograms$sample_ids <- c('ko15', 'wt15')
  opts$chromatograms$aggregation_fun <- 'max'

  # act
  result <- run_create_chromatograms(container, opts, NULL)

  # assert
  expect_gt(nrow(result$chromatograms), 0)
  expect_gte(round(min(result$chromatograms$rt)), 2500)
  expect_lte(round(max(result$chromatograms$rt)), 4500)
  expect_equal(nrow(mass_traces_of(result, opts)), 0)
  expect_equal(nrow(result$detected_peaks), 0)
})

test_that("create_chromatograms with NULL features labels corrected and uncorrected RTs", {
  # arrange
  container <- lcmsPlotDataContainer(
    cached_xcms_object(), sample_id_column = "sample_name")

  opts <- lcmsPlot:::default_options()
  opts$chromatograms$sample_ids <- c('ko15', 'wt15')
  opts$chromatograms$rt_type <- "both"

  # act
  result <- run_create_chromatograms(container, opts, NULL)

  # assert
  expect_equal(
    result$chrom_data$rt_type,
    rep(c("RT uncorrected", "RT corrected"), times = 2))
  by_type <- split(result$chromatograms, result$chromatograms$rt_type)
  expect_equal(
    by_type[["RT corrected"]]$intensity,
    by_type[["RT uncorrected"]]$intensity)
  expect_false(isTRUE(all.equal(
    by_type[["RT corrected"]]$rt, by_type[["RT uncorrected"]]$rt)))
})

test_that("create_chromatograms requires adjusted RTs for corrected chromatograms", {
  # arrange
  container <- lcmsPlotDataContainer(get_test_sample_paths())

  opts <- lcmsPlot:::default_options()
  opts$chromatograms$sample_ids <- c('ko15', 'wt15')
  opts$chromatograms$rt_type <- "corrected"

  # act / assert
  expect_error(
    run_create_chromatograms(container, opts, NULL),
    "adjusted retention times")
})

lipid_search_standards <- function(lipids_query = NULL) {
  # Two standards, but only the kynurenine run is declared in the
  # LipidSearch results; the proline run is an extra sample.
  mzml_dir <- file.path(tempdir(), "lipid-search-standards")
  utils::unzip(
    system.file("extdata", "standards-mzml.zip", package = "lcmsPlot"),
    exdir = mzml_dir
  )
  kynurenine <- file.path(mzml_dir, "mzml", "l-kynurenine-MS1.mzML")
  proline <- file.path(mzml_dir, "mzml", "l-proline-MS1.mzML")

  # Both compounds are reported in the declared sample only. The RT window is
  # widened enough (via rt_extend) to span both runs.
  results_path <- .write_lipid_search_results(
    file_names = "l-kynurenine-MS1.raw",
    lipids = data.frame(
      "Rej." = c(0, 0),
      "LipidIon" = c("KYN(0:0)+H", "PRO(0:0)+H"),
      "LipidGroup" = c("KYN(0:0)", "PRO(0:0)"),
      "Class" = c("KYN", "PRO"),
      "FattyAcid" = c("(0:0)", "(0:0)"),
      "CalcMz" = c(209.0921, 116.0706),
      "IonFormula" = c("C10 H13 O3 N2", "C5 H10 O2 N1"),
      check.names = FALSE
    ),
    values = list(
      "Area" = rbind(24050, 0),
      "Height" = rbind(24050, 0),
      "Rt" = rbind(6.08, 6.08),
      "ObsMz" = rbind("209.0921", ""),
      "Hwhm(L)" = rbind(0.05, 0),
      "Hwhm(R)" = rbind(0.05, 0),
      "Grade" = rbind("A", ""),
      "mScore" = rbind("9", ""),
      "S/N" = rbind("30", "")
    )
  )

  container <- lcmsPlotDataContainer(
    LipidSearchSource(results_path, c(kynurenine, proline)),
    sample_id_column = "sample_id")

  opts <- lcmsPlot:::default_options()
  opts$chromatograms$sample_ids <- sampleMetadata(container)$sample_id
  opts$chromatograms$ppm <- 20
  opts$lipid_search$rt_extend <- 120
  opts$lipid_search$lipids_query <- lipids_query

  list(container = container, opts = opts)
}

test_that("create_chromatograms extracts LipidSearch compounds, including from samples absent from the results", {
  # arrange
  standards <- lipid_search_standards()
  container <- standards$container

  # act
  result <- run_create_chromatograms(container, standards$opts, NULL)

  # assert - both samples yield chromatograms, including the one that the
  # results file never mentioned.
  expect_setequal(unique(result$chromatograms$sample_index), c(1, 2))
  expect_setequal(
    unique(result$chrom_data$name),
    c("KYN(0:0)+H", "PRO(0:0)+H")
  )

  extra_sample_index <- which(!sampleMetadata(container)$in_results)
  extra <- result$chromatograms[
    result$chromatograms$sample_index == extra_sample_index, ]
  expect_gt(nrow(extra), 0)
  expect_gt(max(extra$intensity), 0)

  # Only the declared sample's detected compound is a highlightable peak.
  expect_equal(nrow(result$detected_peaks), 1)
  expect_equal(result$detected_peaks$name, "KYN(0:0)+H")
  expect_false(extra_sample_index %in% result$detected_peaks$sample_index)
})

test_that("create_chromatograms applies the LipidSearch query per lipid, not per lipid/sample row", {
  # arrange - only the kynurenine compound is graded, and only in the one
  # sample the results file declares.
  standards <- lipid_search_standards(lipids_query = 'grade == "A"')
  container <- standards$container
  opts <- standards$opts

  # act
  result <- run_create_chromatograms(container, opts, NULL)

  # assert - the query keeps only the graded lipid, but it is still extracted
  # from every sample, including the one the results file never mentioned.
  expect_equal(unique(result$chrom_data$name), "KYN(0:0)+H")
  expect_setequal(unique(result$chrom_data$sample_index), c(1, 2))

  # A query that matches nothing is an error.
  opts$lipid_search$lipids_query <- 'class == "does-not-exist"'
  expect_error(
    run_create_chromatograms(container, opts, NULL),
    "No lipids matched"
  )
})
