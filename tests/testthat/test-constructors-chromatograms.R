test_that("create_bpc_tic creates a correct base peak chromatogram", {
  # arrange
  raw_data <- get_raw_data()

  # act
  bpc <- create_bpc_tic(raw_data, aggregation_fun = "max")
  close_raw_data(raw_data)

  # assert
  expect_gt(nrow(bpc$chromatograms), 0)

  rt_range <- range(bpc$chromatograms$rt)
  expect_gt(rt_range[1], 2500)
  expect_lt(rt_range[1], 2510)
  expect_gt(rt_range[2], 4400)
  expect_lt(rt_range[2], 4500)
})

test_that("create_bpc_tic supports mean alongside max and sum", {
  # arrange
  raw_data <- get_raw_data()
  hdr <- ms_header(raw_data)
  ms1 <- hdr[hdr$msLevel == 1, ]

  # act
  bpc <- create_bpc_tic(raw_data, aggregation_fun = "max")
  tic <- create_bpc_tic(raw_data, aggregation_fun = "sum")
  aic <- create_bpc_tic(raw_data, aggregation_fun = "mean")
  close_raw_data(raw_data)

  # assert
  expect_equal(bpc$chromatograms$intensity, ms1$basePeakIntensity)
  expect_equal(tic$chromatograms$intensity, ms1$totIonCurrent)
  expect_equal(aic$chromatograms$intensity, ms1$totIonCurrent / ms1$peaksCount)
})

test_that("create_bpc_tic rejects an unknown aggregation function", {
  # arrange
  raw_data <- get_raw_data()

  # act / assert
  expect_error(
    create_bpc_tic(raw_data, aggregation_fun = "median"),
    "Unknown aggregation_fun")
  close_raw_data(raw_data)
})

test_that("ms_header reports a peak count and a mean intensity", {
  # arrange
  raw_data <- get_raw_data()

  # act
  hdr <- ms_header(raw_data)
  close_raw_data(raw_data)

  # assert
  expect_contains(colnames(hdr), c("peaksCount", "meanIntensity"))
  expect_true(all(hdr$meanIntensity >= 0))
})

test_that("lp_chromatogram validates the aggregation function", {
  # act / assert
  expect_error(lp_chromatogram(aggregation_fun = "median"), "should be one of")
})

test_that("create_chromatogram creates a chromatogram within the specified ranges", {
  # arrange
  raw_data <- get_raw_data()

  # act
  chrom <- create_chromatogram(raw_data, mz_range = c(200, 300), rt_range = c(4200, 4300))
  close_raw_data(raw_data)

  # assert
  expect_gt(nrow(chrom$chromatograms), 0)
  expect_gt(nrow(chrom$mass_traces), 0)

  rt_range <- range(chrom$chromatograms$rt)
  expect_gt(rt_range[1], 4200)
  expect_lt(rt_range[2], 4300)

  mz_range <- range(chrom$mass_traces$mz)
  expect_gt(mz_range[1], 200)
  expect_lt(mz_range[2], 300)
})
