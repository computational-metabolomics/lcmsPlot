# TODO: add validation to all slots

#' Create an lcmsPlotDataContainer object from a data object (e.g., XCMSnExp)
#'
#' @param data_obj ...
#' @param sample_id_column ...
#' @param metadata ...
create_data_container_from_obj <- function(data_obj, sample_id_column, metadata) {
  new("lcmsPlotDataContainer",
      data_obj = data_obj,
      metadata = get_metadata(data_obj, sample_id_column, metadata),
      chromatograms = data.frame(),
      mass_traces = data.frame(),
      spectra = data.frame(),
      processed_data_info = data.frame(),
      detected_peaks = data.frame())
}

#' lcmsPlotDataContainer class
#'
#' @slot rts The chromatogram retention time points
#' @slot mzs The mass trace masses
#' @slot intensities The chromatogram intensities
#' @slot sample_ids The sample IDs
#' @slot metadata The metadata for each chromatogram point
#' @export
setClass(
  "lcmsPlotDataContainer",
  slots = list(
    data_obj = "ANY",
    metadata = "data.frame",
    chromatograms = "data.frame",
    mass_traces = "data.frame",
    spectra = "data.frame",
    processed_data_info = "data.frame",
    detected_peaks = "data.frame"
  ),
  prototype = list(
    data_obj = NULL,
    metadata = NULL,
    chromatograms = NULL,
    mass_traces = NULL,
    spectra = NULL,
    processed_data_info = NULL,
    detected_peaks = NULL
  )
)

.create_bpc_tic <- function(raw_data, aggregation_fun) {
  hdr <- mzR::header(raw_data)
  ms1_header <- hdr[hdr$msLevel == 1, ]
  rt <- ms1_header$retentionTime

  if (aggregation_fun == "max") {
    bpi <- ms1_header$basePeakIntensity
  } else {
    bpi <- ms1_header$totIonCurrent
  }

  return(list(
    chromatograms = data.frame(rt = rt, intensity = bpi)
  ))
}

.create_chromatogram <- function(raw_data, mz_range, rt_range) {
  hdr <- mzR::header(raw_data)
  scans_in_rt <- hdr[hdr$retentionTime >= rt_range[1] & hdr$retentionTime <= rt_range[2], ]

  chr <- data.frame()
  mass_traces <- data.frame()

  for (scan_id in scans_in_rt$seqNum) {
    spectrum <- mzR::peaks(raw_data, scan_id)
    rt <- scans_in_rt$retentionTime[scans_in_rt$seqNum == scan_id]
    in_mz_range <- spectrum[spectrum[, 1] >= mz_range[1] & spectrum[, 1] <= mz_range[2], , drop = FALSE]
    total_intensity <- sum(in_mz_range[, 2])

    if (nrow(in_mz_range) > 0) {
      chr <- rbind(chr, data.frame(rt = rt, intensity = total_intensity))
      mass_traces <- rbind(mass_traces, data.frame(rt = rt, mz = in_mz_range[, 1]))
    }
  }

  return(list(
    chromatograms = chr,
    mass_traces = mass_traces
  ))
}

.create_spectrum_from_closest_scan_to_rt <- function(raw_data, rt, ms_level) {
  hdr <- mzR::header(raw_data)

  ms_level_indices <- which(hdr$msLevel == ms_level)

  if (length(ms_level_indices) == 0) {
    # TODO: rethink this
    stop("No scans found with the specified MS level.")
  }

  rt_diffs <- abs(hdr$retentionTime[ms_level_indices] - rt)

  closest_index <- ms_level_indices[which.min(rt_diffs)]

  spectrum_data <- mzR::peaks(raw_data, closest_index)

  spectrum_df <- data.frame(
    mz = spectrum_data[, 1],
    intensity = spectrum_data[, 2],
    rt = hdr$retentionTime[closest_index]
  )

  return(spectrum_df)
}

#' Creates a lcmsPlotDataContainer from a features matrix
#'
#' @param obj A lcmsPlotDataContainer object
#' @param feature_ids ...
#' @param sample_ids ...
#' @param ppm The ppm error for the chromatograms
#' @param rt_tol The RT tolerance for the chromatograms
#' @returns A lcmsPlotDataContainer object
#' @export
setGeneric(
  "create_chromatograms_from_features",
  function(obj, feature_ids, sample_ids, ppm, rt_tol) standardGeneric("create_chromatograms_from_features")
)

#' @rdname create_chromatograms_from_features
setMethod(
  f = "create_chromatograms_from_features",
  signature = c("lcmsPlotDataContainer", "character", "character", "numeric", "numeric"),
  definition = function(obj, feature_ids, sample_ids, ppm, rt_tol) {
    metadata <- obj@metadata %>% filter(sample_id %in% sample_ids)
    raw_data <- io_utils$get_raw_data(metadata$sample_path)
    all_detected_peaks <- get_detected_peaks(obj@data_obj)
    grouped_peaks <- get_grouped_peaks(obj@data_obj) %>% filter(name %in% feature_ids)
    detected_peaks <- data.frame()

    chromatograms <- data.frame()
    mass_traces <- data.frame()
    processed_data_info <- data.frame()

    for (i in seq_len(nrow(grouped_peaks))) {
      feature <- grouped_peaks[i,]
      rtr <- c(feature$rt - rt_tol, feature$rt + rt_tol)
      peak_indices <- as.numeric(unlist(strsplit(feature %>% pull(peakidx), ',')))
      peaks <- all_detected_peaks %>%
        filter(
          row_number() %in% peak_indices,
          sample %in% metadata$sample_index
        ) %>%
        left_join(metadata, by = c("sample" = "sample_index")) # TODO: check this, too specific

      detected_peaks <- rbind(detected_peaks, peaks)

      for (j in seq_len(nrow(peaks))) {
        peak <- peaks[j,]
        mzr <- get_mz_range(peak$mz, ppm)
        sample_metadata <- metadata %>% filter(sample_index == peak$sample)
        raw_obj <- raw_data[[sample_metadata$sample_path]]

        data <- .create_chromatogram(raw_obj, mz_range = mzr, rt_range = rtr)

        processed_data_info <- rbind(processed_data_info, cbind(
          sample_metadata,
          data.frame(
            feature_id = feature$name,
            peak_rt_min = peak$rtmin,
            peak_rt_max = peak$rtmax
          )
        ))

        chromatograms <- rbind(chromatograms, data.frame(
          rt = data$chromatograms$rt,
          intensity = data$chromatograms$intensity,
          metadata_index = nrow(processed_data_info)
        ))

        mass_traces <- rbind(mass_traces, data.frame(
          rt = data$mass_traces$rt,
          mz = data$mass_traces$mz,
          metadata_index = nrow(processed_data_info)
        ))
      }
    }

    obj@chromatograms <- chromatograms
    obj@mass_traces <- mass_traces
    obj@processed_data_info <- processed_data_info
    obj@detected_peaks = detected_peaks

    return(obj)
  }
)

#' Creates a lcmsPlotDataContainer from raw sample files
#'
#' @param raw_data A list of mzR objects
#' @param features A matrix of entries with mz and rt values
#' @param sample_ids The sample IDs to select
#' @param ppm The ppm error for the chromatograms
#' @param rt_tol The RT tolerance for the chromatograms
#' @returns A lcmsPlotDataContainer object
#' @export
setGeneric(
  "create_chromatograms_from_raw",
  function(obj, features, sample_ids, ppm, rt_tol) standardGeneric("create_chromatograms_from_raw")
)

#' @rdname create_chromatograms_from_raw
setMethod(
  f = "create_chromatograms_from_raw",
  signature = c("lcmsPlotDataContainer", "matrix", "character", "numeric", "numeric"),
  definition = function(obj, features, sample_ids, ppm, rt_tol) {
    metadata <- obj@metadata %>% filter(sample_id %in% sample_ids)
    raw_data <- io_utils$get_raw_data(metadata$sample_path)
    all_detected_peaks <- get_detected_peaks(obj@data_obj)

    chromatograms <- data.frame()
    mass_traces <- data.frame()
    processed_data_info <- data.frame()
    detected_peaks <- data.frame()

    for (i in 1:nrow(features)) {
      feature <- features[i,]
      feature_id <- paste0('M', round(feature['mz']), 'T', round(feature['rt']))
      mzr <- get_mz_range(feature['mz'], ppm)
      rtr <- c(feature['rt'] - rt_tol, feature['rt'] + rt_tol)

      for (j in 1:nrow(metadata)) {
        sample_metadata <- metadata[j,]
        raw_obj <- raw_data[[sample_metadata$sample_path]]

        data <- .create_chromatogram(raw_obj, mz_range = mzr, rt_range = rtr)

        if (!is.null(all_detected_peaks)) {
          peaks <- all_detected_peaks %>%
            filter(
              sample %in% sample_metadata$sample_index,
              mz >= mzr[1],
              mz <= mzr[2],
              rt >= rtr[1],
              rt <= rtr[2]
            ) %>%
            left_join(metadata, by = c("sample" = "sample_index")) # TODO: check this, too specific
        } else {
          peaks <- data.frame()
        }

        detected_peaks <- rbind(detected_peaks, peaks)

        processed_data_info <- rbind(processed_data_info, cbind(
          sample_metadata,
          data.frame(feature_id = feature_id)
        ))

        chromatograms <- rbind(chromatograms, data.frame(
          rt = data$chromatograms$rt,
          intensity = data$chromatograms$intensity,
          metadata_index = nrow(processed_data_info)
        ))

        mass_traces <- rbind(mass_traces, data.frame(
          rt = data$mass_traces$rt,
          mz = data$mass_traces$mz,
          metadata_index = nrow(processed_data_info)
        ))
      }
    }

    obj@chromatograms <- chromatograms
    obj@mass_traces <- mass_traces
    obj@processed_data_info <- processed_data_info
    obj@detected_peaks <- detected_peaks

    return(obj)
  }
)

#' Creates a lcmsPlotDataContainer from raw sample files
#'
#' @param obj A lcmsPlotDataContainer object
#' @param sample_ids The sample IDs to select
#' @param aggregation_fun Specify the function to be used to aggregate intensity
#' values across the mz value range for the same retention time.
#' Allowed values are "sum" (the default), "max"
#' @returns A lcmsPlotDataContainer object
#' @export
setGeneric(
  "create_full_rt_chromatograms",
  function(obj, sample_ids, aggregation_fun) standardGeneric("create_full_rt_chromatograms")
)

#' @rdname create_full_rt_chromatograms
setMethod(
  f = "create_full_rt_chromatograms",
  signature = c("lcmsPlotDataContainer", "character", "character"),
  definition = function(obj, sample_ids, aggregation_fun) {
    metadata <- obj@metadata %>% filter(sample_id %in% sample_ids)
    raw_data <- io_utils$get_raw_data(metadata$sample_path)

    chromatograms <- data.frame()
    processed_data_info <- data.frame()

    for (i in 1:nrow(metadata)) {
      sample_metadata <- metadata[i,]
      raw_obj <- raw_data[[sample_metadata$sample_path]]

      data <- .create_bpc_tic(raw_obj, aggregation_fun)

      processed_data_info <- rbind(processed_data_info, sample_metadata)

      chromatograms <- rbind(chromatograms, data.frame(
        rt = data$chromatograms$rt,
        intensity = data$chromatograms$intensity,
        metadata_index = nrow(processed_data_info)
      ))
    }

    obj@chromatograms <- chromatograms
    obj@processed_data_info <- processed_data_info

    return(obj)
  }
)



#' Creates/updates am lcmsPlotDataContainer from spectra
#'
#' @param obj A lcmsPlotDataContainer object
#' @param options ...
#' @returns A lcmsPlotDataContainer object
#' @export
setGeneric(
  "create_spectra",
  function(obj, options) standardGeneric("create_spectra")
)

#' @rdname create_spectra
setMethod(
  f = "create_spectra",
  signature = c("lcmsPlotDataContainer", "list"),
  definition = function(obj, options) {
    metadata <- obj@metadata %>% filter(sample_id %in% options$sample_ids)

    all_spectra <- data.frame()

    for (i in seq_len(nrow(metadata))) {
      sample_metadata <- metadata[i,]
      raw_obj <- mzR::openMSfile(sample_metadata$sample_path)

      if (options$spectra$mode == "closest" & !is.null(options$spectra$rt)) {
        spectra <- .create_spectrum_from_closest_scan_to_rt(raw_obj, options$spectra$rt, options$spectra$ms_level)
      } else if (options$spectra$mode == "closest_apex") {
        spectra <- obj@detected_peaks %>%
          filter(sample_id == sample_metadata$sample_id) %>%
          pull(rt) %>%
          lapply(function (rt) .create_spectrum_from_closest_scan_to_rt(raw_obj, rt, options$spectra$ms_level)) %>%
          do.call(rbind, .)
      } else if (options$spectra$mode == "across_peak") {
        spectra <- expanded_peaks <- obj@detected_peaks %>%
          rowwise() %>%
          mutate(intervals = list(seq(rtmin, rtmax, by = options$spectra$interval))) %>%
          tidyr::unnest(intervals) %>%
          mutate(rt_interval = intervals) %>%
          select(sample_id, rt, rtmin, rtmax, rt_interval) %>%
          pull(rt_interval) %>%
          lapply(function (rt) .create_spectrum_from_closest_scan_to_rt(raw_obj, rt, options$spectra$ms_level)) %>%
          do.call(rbind, .)

      } else {
        # TODO: throw an error
      }

      spectra <- spectra %>%
        mutate(metadata_index = i) # TODO: check processed_data_info
      all_spectra <- rbind(all_spectra, spectra)

      mzR::close(raw_obj)
    }

    obj@spectra <- all_spectra

    return(obj)
  }
)
