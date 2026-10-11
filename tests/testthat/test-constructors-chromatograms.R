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

# A fake rawrr backend: XICs are a Gaussian centred on each mass's integer
# part (in minutes), and every call is counted.
local_fake_rawrr <- function(env = parent.frame()) {
  calls <- new.env()
  calls$xic <- list()
  calls$index <- 0
  calls$spectrum <- 0
  times <- seq(0, 20, by = 0.05)

  local_mocked_bindings(
    readChromatogram = function(rawfile, mass = NULL, tol = 10, type = "xic", ...) {
      if (type == "xic") {
        calls$xic[[length(calls$xic) + 1]] <- list(mass = mass, tol = tol)
        return(lapply(mass, function(m) {
          list(times = times, intensities = 1e5 * dnorm(times, m %% 20, 0.2))
        }))
      }
      list(times = times, intensities = rep(1, length(times)))
    },
    readIndex = function(rawfile) {
      calls$index <- calls$index + 1
      data.frame(scan = seq_along(times), StartTime = times, MSOrder = "Ms")
    },
    readSpectrum = function(rawfile, scan) {
      calls$spectrum <- calls$spectrum + 1
      lapply(scan, function(s) list(mZ = c(100, 201, 300), intensity = c(1, 2, 3)))
    },
    .package = "rawrr",
    .env = env
  )

  calls
}

new_fake_rawrr_reader <- function() {
  new("RawrrReader", path = "fake.raw", cache = new.env(parent = emptyenv()))
}

test_that("create_chromatograms_batch reads all rawrr windows in one call", {
  skip_if_not_installed("rawrr")
  calls <- local_fake_rawrr()
  reader <- new_fake_rawrr_reader()
  mzs <- c(201, 305, 412)
  mz_ranges <- do.call(rbind, lapply(mzs, get_mz_range, ppm = 5))
  rt_ranges <- rbind(c(30, 90), c(270, 330), c(700, 760))

  # act
  batch <- create_chromatograms_batch(reader, mz_ranges, rt_ranges)

  # assert
  expect_length(calls$xic, 1)
  expect_equal(calls$xic[[1]]$mass, mzs)
  expect_length(batch, 3)
  for (i in seq_along(mzs)) {
    single <- create_chromatogram(
      reader, mz_ranges[i, ], rt_ranges[i, ], include_mass_traces = FALSE)
    expect_equal(batch[[i]], single$chromatograms)
    expect_true(all(batch[[i]]$rt >= rt_ranges[i, 1]))
    expect_true(all(batch[[i]]$rt <= rt_ranges[i, 2]))
  }
})

test_that("ms_chromatogram makes one rawrr call per distinct tolerance", {
  skip_if_not_installed("rawrr")
  calls <- local_fake_rawrr()
  reader <- new_fake_rawrr_reader()

  # act
  chroms <- ms_chromatogram(
    reader, mz = c(201, 305, 412), ppm = c(10, 20, 10), rt_range = c(0, 1200))

  # assert
  expect_length(calls$xic, 2)
  expect_setequal(
    vapply(calls$xic, `[[`, numeric(1), "tol"), c(10, 20))
  peak_rt <- vapply(chroms, function(x) x$rt[which.max(x$intensity)], numeric(1))
  expect_equal(peak_rt, c(1, 5, 12) * 60)
})

test_that("create_chromatograms_batch returns empty chromatograms for NA windows", {
  skip_if_not_installed("rawrr")
  calls <- local_fake_rawrr()
  reader <- new_fake_rawrr_reader()

  # act
  batch <- create_chromatograms_batch(
    reader,
    mz_ranges = rbind(get_mz_range(201, 5), c(NA, NA)),
    rt_ranges = rbind(c(30, 90), c(30, 90)))

  # assert
  expect_equal(calls$xic[[1]]$mass, 201)
  expect_gt(nrow(batch[[1]]), 0)
  expect_equal(nrow(batch[[2]]), 0)
})

test_that("ms_header is read once per rawrr reader", {
  skip_if_not_installed("rawrr")
  calls <- local_fake_rawrr()
  reader <- new_fake_rawrr_reader()

  # act
  first <- ms_header(reader)
  second <- ms_header(reader)

  # assert
  expect_equal(calls$index, 1)
  expect_identical(first, second)
})

test_that("create_chromatogram skips rawrr header and spectra without mass traces", {
  skip_if_not_installed("rawrr")
  calls <- local_fake_rawrr()
  reader <- new_fake_rawrr_reader()
  mzr <- get_mz_range(201, 5)

  # act
  without <- create_chromatogram(reader, mzr, c(30, 90), include_mass_traces = FALSE)

  # assert
  expect_equal(calls$index, 0)
  expect_equal(calls$spectrum, 0)
  expect_equal(nrow(without$mass_traces), 0)

  # act
  with <- create_chromatogram(reader, mzr, c(30, 90))

  # assert
  expect_equal(calls$index, 1)
  expect_equal(calls$spectrum, 1)
  expect_gt(nrow(with$mass_traces), 0)
  expect_equal(with$chromatograms, without$chromatograms)
})

test_that("create_chromatograms_batch matches create_chromatogram on mzR files", {
  # arrange
  raw_data <- get_raw_data()
  mz_ranges <- rbind(c(200, 300), c(300.1, 300.4), c(NA, NA))
  rt_ranges <- rbind(c(4200, 4300), c(3000, 3200), c(3000, 3200))

  # act
  batch <- create_chromatograms_batch(
    raw_data, mz_ranges, rt_ranges, fill_gaps = TRUE)
  singles <- lapply(1:2, function(i) {
    create_chromatogram(
      raw_data, mz_ranges[i, ], rt_ranges[i, ], fill_gaps = TRUE
    )$chromatograms
  })
  close_raw_data(raw_data)

  # assert
  expect_equal(batch[1:2], singles)
  expect_equal(nrow(batch[[3]]), 0)
})

test_that("ms_header is memoised per mzR reader", {
  # arrange
  raw_data <- get_raw_data()

  # act
  first <- ms_header(raw_data)
  second <- ms_header(raw_data)

  # assert
  expect_identical(first, second)
  expect_identical(raw_data@cache$header, first)
  close_raw_data(raw_data)
})

test_that("mzR readers do not share a cache", {
  # arrange
  reader_a <- get_raw_data(1)
  reader_b <- get_raw_data(2)

  # act
  ms_header(reader_a)

  # assert
  expect_false(identical(reader_a@cache, reader_b@cache))
  expect_null(reader_b@cache$header)
  close_raw_data(reader_a)
  close_raw_data(reader_b)
})

test_that("ms_scan_index matches the header without reading the spectra", {
  # arrange
  raw_data <- get_raw_data()

  # act
  scan_index <- ms_scan_index(raw_data)

  # assert: on netCDF the full header is what reads every spectrum
  expect_null(raw_data@cache$header)
  hdr <- ms_header(raw_data)
  expect_equal(
    as.data.frame(scan_index),
    as.data.frame(hdr[, c("seqNum", "retentionTime", "msLevel")]))
  close_raw_data(raw_data)
})

# The scan-by-scan definition of an extracted ion chromatogram, kept as a
# reference for the vectorised implementation.
naive_chromatogram <- function(raw_data, mz_range, rt_range, fill_gaps) {
  hdr <- ms_header(raw_data)
  hdr <- hdr[hdr$msLevel == 1, ]
  hdr <- hdr[hdr$retentionTime >= rt_range[1] &
               hdr$retentionTime <= rt_range[2], ]
  chr <- data.frame(rt = numeric(), intensity = numeric())
  mt <- data.frame(rt = numeric(), mz = numeric())
  for (i in seq_len(nrow(hdr))) {
    pk <- ms_peaks(raw_data, hdr$seqNum[i])[[1]]
    pk <- pk[pk[, 1] >= mz_range[1] & pk[, 1] <= mz_range[2], , drop = FALSE]
    if (nrow(pk) > 0) {
      chr <- rbind(chr, data.frame(rt = hdr$retentionTime[i], intensity = sum(pk[, 2])))
      mt <- rbind(mt, data.frame(rt = hdr$retentionTime[i], mz = pk[, 1]))
    } else if (fill_gaps) {
      chr <- rbind(chr, data.frame(rt = hdr$retentionTime[i], intensity = 0))
    }
  }
  rownames(chr) <- NULL
  rownames(mt) <- NULL
  list(chromatograms = chr, mass_traces = mt)
}

test_that("create_chromatogram matches the scan-by-scan reference", {
  # arrange
  raw_data <- get_raw_data()
  cases <- list(
    list(mz = c(334.9, 335.1), rt = c(2700, 2900), fill_gaps = FALSE),
    list(mz = c(334.9, 335.1), rt = c(2700, 2900), fill_gaps = TRUE),
    list(mz = c(200, 300), rt = c(4200, 4300), fill_gaps = FALSE))

  for (case in cases) {
    # act
    actual <- create_chromatogram(
      raw_data, case$mz, case$rt, fill_gaps = case$fill_gaps)
    expected <- naive_chromatogram(raw_data, case$mz, case$rt, case$fill_gaps)

    # assert
    expect_equal(as.data.frame(actual$chromatograms), expected$chromatograms)
    expect_equal(as.data.frame(actual$mass_traces), expected$mass_traces)
  }
  close_raw_data(raw_data)
})

test_that("create_chromatogram keeps its columns for an empty window", {
  # arrange
  raw_data <- get_raw_data()

  # act: an RT window outside the run
  outside <- create_chromatogram(raw_data, c(334.9, 335.1), c(1, 2))
  # act: scans in range but no m/z in the window, with gap filling
  no_mz <- create_chromatogram(
    raw_data, c(5000, 5001), c(2700, 2710), fill_gaps = TRUE)
  close_raw_data(raw_data)

  # assert
  expect_equal(nrow(outside$chromatograms), 0)
  expect_named(outside$mass_traces, c("rt", "mz"))
  expect_equal(nrow(outside$mass_traces), 0)
  expect_gt(nrow(no_mz$chromatograms), 0)
  expect_true(all(no_mz$chromatograms$intensity == 0))
  expect_named(no_mz$mass_traces, c("rt", "mz"))
  expect_equal(nrow(no_mz$mass_traces), 0)
})
