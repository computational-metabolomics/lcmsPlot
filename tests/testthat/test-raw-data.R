test_that("XcmsRawList data is converted into per-scan spectra", {
  # arrange
  raw <- xcms::xcmsRaw(get_test_sample_paths(1), profstep = 0)

  # act
  sps <- xcmsraw_list_to_spectra(XcmsRawList(raw))

  # assert
  expect_length(sps, length(raw@scantime))
  expect_equal(Spectra::rtime(sps), raw@scantime)
  expect_equal(sum(lengths(sps)), length(raw@env$mz))
  scan_10 <- (raw@scanindex[10] + 1):raw@scanindex[11]
  expect_equal(Spectra::mz(sps)[[10]], raw@env$mz[scan_10])
  expect_equal(Spectra::intensity(sps)[[10]], raw@env$intensity[scan_10])
  expect_equal(sps$sample_index, rep(1L, length(sps)))
})
