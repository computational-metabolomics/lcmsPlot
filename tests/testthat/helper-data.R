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
    pdata = new("NAnnotatedDataFrame", pd),
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
