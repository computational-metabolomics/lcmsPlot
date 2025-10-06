get_raw_data <- function(file_index = 1) {
  raw_file <- dir(system.file("cdf", package = "faahKO"), full.names = TRUE, recursive = TRUE)[file_index]
  mzR::openMSfile(raw_file)
}

close_raw_data <- function(raw_data) {
  mzR::close(raw_data)
}

get_test_sample_paths <- function(indices = c(1, 7)) {
  dir(system.file("cdf", package = "faahKO"), full.names = TRUE, recursive = TRUE)[indices]
}

get_MsExperiment_object <- function() {
  cdfs <- dir(system.file("cdf", package = "faahKO"), full.names = TRUE,
              recursive = TRUE)[c(1, 7)]
  sample_names <- sub(basename(cdfs), pattern = ".CDF", replacement = "", fixed = TRUE)

  pd <- data.frame(
    sample_name = sample_names,
    sample_group = c("KO", "WT"),
    stringsAsFactors = FALSE)

  MsExperiment::readMsExperiment(spectraFiles = cdfs, sampleData = pd)
}

get_XCMSnExp_object <- function(should_detect_peaks = FALSE, should_group_peaks = FALSE) {
  cdfs <- dir(system.file("cdf", package = "faahKO"), full.names = TRUE,
              recursive = TRUE)[c(1, 7)]
  sample_names <- sub(basename(cdfs), pattern = ".CDF", replacement = "", fixed = TRUE)

  pd <- data.frame(
    sample_name = sample_names,
    sample_group = c("KO", "WT"),
    stringsAsFactors = FALSE)

  raw_data <- MSnbase::readMSData(
    files = cdfs,
    pdata = new("AnnotatedDataFrame", pd),
    mode = "onDisk",
    msLevel = 1)

  xdata <- as(raw_data, "XCMSnExp")

  if (should_detect_peaks) {
    cwp <- CentWaveParam(peakwidth = c(20, 80), noise = 10000, prefilter = c(6, 10000))
    xdata <- findChromPeaks(xdata, param = cwp)
  }

  if (should_detect_peaks && should_group_peaks) {
    xdata <- adjustRtime(xdata, param = ObiwarpParam(binSize = 0.6))
    pdp <- PeakDensityParam(
      sampleGroups = pd$sample_group,
      minFraction = 1,
      bw = 30)
    xdata <- groupChromPeaks(xdata, param = pdp)
  }

  xdata
}

generate_chromatogram <- function(sample_id) {
  rt <- seq(0, 20, by = 0.1)

  peak1 <- dnorm(rt, mean = 5, sd = 0.3) * 100
  peak2 <- dnorm(rt, mean = 12, sd = 0.5) * 150
  peak3 <- dnorm(rt, mean = 17, sd = 0.4) * 120

  noise <- runif(length(rt), min = 5, max = 15)

  intensity <- peak1 + peak2 + peak3 + noise

  data.frame(rt = rt, intensity = intensity, sample_id = sample_id)
}

generate_spectra <- function(sample_id) {
  mz <- sort(runif(100, min = 50, max = 1000))
  intensity <- abs(rnorm(100, mean = 1e5, sd = 5e4))
  data.frame(
    mz = mz,
    intensity = intensity,
    rt = runif(1, min = 50, max = 1000),
    reference = FALSE,
    sample_id = sample_id
  )
}

generate_datasets_for_plots <- function(n_samples = 5) {
  chromatograms <- data.frame()
  spectra <- data.frame()

  for (i in seq_len(n_samples)) {
    chromatograms <- rbind(
      chromatograms,
      generate_chromatogram(paste0("sample_", i))
    )

    spectra <- rbind(
      spectra,
      generate_spectra(paste0("sample_", i))
    )
  }

  return(list(
    chromatograms = chromatograms,
    spectra = spectra
  ))
}

get_facets_from_plot <- function(plot) {
  f <- plot$facet
  if (inherits(f, "FacetNull")) {
    return(NULL)  # no facets
  } else if (inherits(f, "FacetWrap")) {
    return(as.character(f$facets))  # names of facetting variables
  } else if (inherits(f, "FacetGrid")) {
    return(c(names(f$rows), names(f$cols)))
  } else {
    return(NULL)  # fallback
  }
}
