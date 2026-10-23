# The scan-by-scan definitions of the chromatograms, kept as a reference for the
# batched implementation. They read the peaks of one sample's MS1 scans.
reference_scans <- function(container, sample_index = 1) {
  sps <- raw_spectra(container, sample_index, ms_level = 1L)
  list(
    rt = Spectra::rtime(sps),
    peaks = as.list(Spectra::peaksData(sps)))
}

naive_chromatogram <- function(scans, mz_range, rt_range, fill_gaps) {
  chr <- data.frame(rt = numeric(), intensity = numeric())
  mt <- data.frame(rt = numeric(), mz = numeric())
  for (i in which(scans$rt >= rt_range[1] & scans$rt <= rt_range[2])) {
    pk <- scans$peaks[[i]]
    pk <- pk[pk[, 1] >= mz_range[1] & pk[, 1] <= mz_range[2], , drop = FALSE]
    if (nrow(pk) > 0) {
      chr <- rbind(chr, data.frame(rt = scans$rt[i], intensity = sum(pk[, 2])))
      mt <- rbind(mt, data.frame(rt = scans$rt[i], mz = pk[, 1]))
    } else if (fill_gaps) {
      chr <- rbind(chr, data.frame(rt = scans$rt[i], intensity = 0))
    }
  }
  rownames(chr) <- NULL
  rownames(mt) <- NULL
  list(chromatograms = chr, mass_traces = mt)
}

one_window <- function(mz_range, rt_range, sample_index = 1L) {
  tibble::tibble(
    sample_index = sample_index,
    feature_id = "F1",
    mzMin = mz_range[1],
    mzMax = mz_range[2],
    rtMin = rt_range[1],
    rtMax = rt_range[2])
}

test_that("create_total_chromatograms creates a base peak chromatogram over the run", {
  # arrange
  container <- get_raw_container()

  # act
  container <- create_total_chromatograms(
    container, .sample_table(container), aggregation_fun = "max")
  bpc <- as_plot_df(container@chromatograms)

  # assert
  expect_gt(nrow(bpc), 0)
  rt_range <- range(bpc$rt)
  expect_gt(rt_range[1], 2500)
  expect_lt(rt_range[1], 2510)
  expect_gt(rt_range[2], 4400)
  expect_lt(rt_range[2], 4500)
})

test_that("create_total_chromatograms summarises each scan's peaks", {
  # arrange
  container <- get_raw_container()
  scans <- reference_scans(container)
  samples <- .sample_table(container)

  # act
  intensity <- function(fun) {
    chromatograms <- create_total_chromatograms(
      container, samples, aggregation_fun = fun)@chromatograms
    as_plot_df(chromatograms)$intensity
  }

  # assert: netCDF headers are unreliable, so every value comes from the peaks
  expect_equal(
    intensity("max"),
    vapply(scans$peaks, function(p) max(p[, 2]), numeric(1)))
  expect_equal(
    intensity("sum"),
    vapply(scans$peaks, function(p) sum(p[, 2]), numeric(1)))
  expect_equal(
    intensity("mean"),
    vapply(scans$peaks, function(p) mean(p[, 2]), numeric(1)))
})

test_that("create_total_chromatograms rejects an unknown aggregation function", {
  # arrange
  container <- get_raw_container()

  # act / assert
  expect_error(
    create_total_chromatograms(
      container, .sample_table(container), aggregation_fun = "median"),
    "should be one of")
})

test_that("create_total_chromatograms takes base peaks and TIC from mzML headers", {
  # arrange: an mzML header reports both values, so no peak is read
  sps <- make_spectra(
    S4Vectors::DataFrame(
      msLevel = 1L,
      rtime = c(1, 2),
      dataOrigin = "a.mzML",
      basePeakIntensity = c(10, 20),
      totIonCurrent = c(100, 200)),
    mz = list(c(1, 2), c(3, 4)),
    intensity = list(c(1, 1), c(1, 1)))

  # act / assert
  expect_equal(lcmsPlot:::.scan_intensity(sps, "max"), c(10, 20))
  expect_equal(lcmsPlot:::.scan_intensity(sps, "sum"), c(100, 200))
  expect_equal(lcmsPlot:::.scan_intensity(sps, "mean"), c(50, 100))
})

test_that("lp_chromatogram validates the aggregation function", {
  # act / assert
  expect_error(lp_chromatogram(aggregation_fun = "median"), "should be one of")
})

test_that("extract_chromatograms matches the scan-by-scan reference", {
  # arrange
  container <- get_raw_container()
  scans <- reference_scans(container)
  cases <- list(
    list(mz = c(334.9, 335.1), rt = c(2700, 2900), fill_gaps = FALSE),
    list(mz = c(334.9, 335.1), rt = c(2700, 2900), fill_gaps = TRUE),
    list(mz = c(200, 300), rt = c(4200, 4300), fill_gaps = FALSE))

  for (case in cases) {
    window <- one_window(case$mz, case$rt)

    # act
    chromatograms <- extract_chromatograms(
      container, window, fill_gaps = case$fill_gaps)@chromatograms
    points <- extract_points(container, window)
    expected <- naive_chromatogram(scans, case$mz, case$rt, case$fill_gaps)

    # assert
    actual <- as_plot_df(chromatograms)
    expect_equal(
      as.data.frame(actual[, c("rt", "intensity")]),
      expected$chromatograms)
    expect_equal(actual$feature_id, rep("F1", nrow(actual)))
    expect_equal(
      as.data.frame(points[, c("rt", "mz")]),
      expected$mass_traces)
  }
})

test_that("extract_chromatograms reads all windows of a sample in one pass", {
  # arrange
  container <- get_raw_container()
  scans <- reference_scans(container)
  windows <- dplyr::bind_rows(
    one_window(c(334.9, 335.1), c(2700, 2900)),
    one_window(c(200, 300), c(4200, 4300)),
    one_window(c(NA, NA), c(2700, 2900)))
  windows$feature_id <- c("A", "B", "C")

  # act
  chromatograms <- extract_chromatograms(container, windows)@chromatograms

  # assert
  expect_length(chromatograms, 3)
  expect_equal(Chromatograms::chromData(chromatograms)$feature_id, c("A", "B", "C"))
  peaks <- Chromatograms::peaksData(chromatograms)
  expect_equal(
    peaks[[2]]$intensity,
    naive_chromatogram(scans, c(200, 300), c(4200, 4300), FALSE)$chromatograms$intensity)
  expect_equal(nrow(peaks[[3]]), 0)
})

test_that("extract_chromatograms keeps an empty window empty", {
  # arrange
  container <- get_raw_container()

  # act: an RT window outside the run
  outside <- extract_chromatograms(
    container, one_window(c(334.9, 335.1), c(1, 2)))@chromatograms
  # act: scans in range but no m/z in the window, with gap filling
  no_mz <- extract_chromatograms(
    container,
    one_window(c(5000, 5001), c(2700, 2710)),
    fill_gaps = TRUE)@chromatograms

  # assert
  expect_equal(nrow(Chromatograms::peaksData(outside)[[1]]), 0)
  no_mz_peaks <- Chromatograms::peaksData(no_mz)[[1]]
  expect_gt(nrow(no_mz_peaks), 0)
  expect_true(all(no_mz_peaks$intensity == 0))
  expect_equal(
    nrow(extract_points(container, one_window(c(334.9, 335.1), c(1, 2)))),
    0)
})

test_that("extract_chromatograms uses the adjusted retention times when asked", {
  # arrange
  container <- lcmsPlotDataContainer(
    cached_xcms_object(), sample_id_column = "sample_name")
  container <- ensure_spectra(container)
  window <- one_window(c(334.9, 335.1), c(2700, 2900), sample_index = 2L)
  sps <- raw_spectra(container, 2L, ms_level = 1L)

  # act
  raw <- as_plot_df(extract_chromatograms(container, window)@chromatograms)
  adjusted <- as_plot_df(extract_chromatograms(
    container, window, adjusted = TRUE)@chromatograms)

  # assert
  expect_true(all(raw$rt %in% Spectra::rtime(sps)))
  expect_true(all(adjusted$rt %in% sps$rtime_adjusted))
  expect_false(isTRUE(all.equal(raw$rt, adjusted$rt)))
})

test_that("extract_chromatograms gives the same result through BiocParallel", {
  # arrange
  container <- get_raw_container(c(1, 2))
  windows <- dplyr::bind_rows(
    one_window(c(334.9, 335.1), c(2700, 2900), sample_index = 1L),
    one_window(c(334.9, 335.1), c(2700, 2900), sample_index = 2L))

  # act
  serial <- extract_chromatograms(container, windows)@chromatograms
  parallel <- extract_chromatograms(
    container, windows, BPPARAM = BiocParallel::SerialParam())@chromatograms

  # assert
  expect_equal(as_plot_df(parallel), as_plot_df(serial))
})

# A fake rawrr backend: XICs are a Gaussian centred on each mass's integer
# part (in minutes), and every call is counted.
local_fake_rawrr <- function(env = parent.frame()) {
  calls <- new.env()
  calls$xic <- list()
  calls$index <- 0
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
    .package = "rawrr",
    .env = env
  )

  calls
}

test_that("rawrr_xic makes one rawrr call per distinct tolerance", {
  skip_if_not_installed("rawrr")
  calls <- local_fake_rawrr()
  mzs <- c(201, 305, 412)
  ppms <- c(10, 20, 10)
  mz_ranges <- do.call(rbind, Map(get_mz_range, mzs, ppms))
  rt_ranges <- matrix(c(0, 1200), nrow = 3, ncol = 2, byrow = TRUE)

  # act
  chroms <- rawrr_xic("fake.raw", mz_ranges, rt_ranges)

  # assert
  expect_length(calls$xic, 2)
  expect_setequal(
    round(vapply(calls$xic, `[[`, numeric(1), "tol")), c(20, 40))
  peak_rt <- vapply(chroms, function(x) x$rtime[which.max(x$intensity)], numeric(1))
  expect_equal(peak_rt, c(1, 5, 12) * 60)
})

test_that("rawrr_xic keeps each window's RT range and skips undefined windows", {
  skip_if_not_installed("rawrr")
  calls <- local_fake_rawrr()
  mz_ranges <- rbind(get_mz_range(201, 5), c(NA, NA), get_mz_range(305, 5))
  rt_ranges <- rbind(c(30, 90), c(30, 90), c(270, 330))

  # act
  chroms <- rawrr_xic("fake.raw", mz_ranges, rt_ranges)

  # assert
  expect_length(calls$xic, 1)
  expect_equal(calls$xic[[1]]$mass, c(201, 305))
  expect_equal(nrow(chroms[[2]]), 0)
  expect_true(all(chroms[[1]]$rtime >= 30 & chroms[[1]]$rtime <= 90))
  expect_true(all(chroms[[3]]$rtime >= 270 & chroms[[3]]$rtime <= 330))
})

test_that("extract_chromatograms reads .raw windows with rawrr, without spectra", {
  skip_if_not_installed("rawrr")
  calls <- local_fake_rawrr()
  container <- lcmsPlotDataContainer(c("a.raw", "b.raw"))
  windows <- dplyr::bind_rows(
    one_window(get_mz_range(201, 5), c(30, 90), sample_index = 1L),
    one_window(get_mz_range(305, 5), c(270, 330), sample_index = 1L),
    one_window(get_mz_range(201, 5), c(30, 90), sample_index = 2L))

  # act
  container <- extract_chromatograms(container, windows)

  # assert: one rawrr call per file, and no spectra were loaded
  expect_length(calls$xic, 2)
  expect_length(container@chromatograms, 3)
  expect_length(MsExperiment::spectra(container@experiment), 0)
})

test_that("create_total_chromatograms reads .raw files with rawrr", {
  skip_if_not_installed("rawrr")
  calls <- local_fake_rawrr()
  container <- lcmsPlotDataContainer("a.raw")

  # act
  bpc <- create_total_chromatograms(
    container, .sample_table(container), aggregation_fun = "max")

  # assert
  expect_equal(calls$index, 1)
  expect_equal(nrow(as_plot_df(bpc@chromatograms)), 401)
  expect_error(
    create_total_chromatograms(
      container, .sample_table(container), aggregation_fun = "mean"),
    "not available")
})
