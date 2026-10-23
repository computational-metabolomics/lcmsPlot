test_that("subset_features keeps the chromatograms and peaks of one feature", {
  # arrange
  p <- lcmsPlot(cached_xcms_object(), sample_id_column = "sample_name") +
    lp_chromatogram(features = baseline_feature_ids(), highlight_peaks = TRUE) +
    lp_spectra(mode = "closest_apex")
  feature <- baseline_feature_ids()[2]

  # act
  one <- subset_features(p, feature_id == feature)

  # assert
  chrom_data <- Chromatograms::chromData(chromatograms(one))
  expect_gt(nrow(chrom_data), 0)
  expect_true(all(chrom_data$feature_id == feature))
  expect_true(all(detectedPeaks(one)$feature_id == feature))
  expect_true(all(spectra(one)$feature_id == feature))
  expect_equal(
    nrow(chrom_data),
    sum(Chromatograms::chromData(chromatograms(p))$feature_id == feature))
})

test_that("subset_features leaves datasets without the columns untouched", {
  # arrange
  p <- lcmsPlot(get_test_sample_paths()) +
    lp_chromatogram(features = baseline_feature_windows()) +
    lp_total_ion_current()
  feature <- Chromatograms::chromData(chromatograms(p))$feature_id[1]

  # act
  one <- subset_features(p, feature_id == feature)

  # assert
  expect_length(chromatograms(one), 2)
  expect_identical(
    summaries(one)$total_ion_current, summaries(p)$total_ion_current)
})
