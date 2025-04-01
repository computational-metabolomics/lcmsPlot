#' Create a LCMSPlotDataContainer object
#'
#' @param chromatograms_df A data frame containing columns intensity, rt, sample_id, and any additional metadata columns
#' @returns A LCMSPlotDataContainer object
#' @export
peakly_chromatogram_container <- function(chromatograms_df, mass_traces_df, metadata_df) {
  new("LCMSPlotDataContainer",
      chromatograms = chromatograms_df,
      mass_traces = mass_traces_df,
      metadata = metadata_df
  )
}

#' LCMSPlotDataContainer class
#'
#' @slot rts The chromatogram retention time points
#' @slot mzs The mass trace masses
#' @slot intensities The chromatogram intensities
#' @slot sample_ids The sample IDs
#' @slot metadata The metadata for each chromatogram point
#' @export
setClass(
  "LCMSPlotDataContainer",
  slots = list(
    chromatograms = "data.frame",
    mass_traces = "data.frame",
    metadata = "data.frame"
  ),
  prototype = list(
    chromatograms = NULL,
    mass_traces = NULL,
    metadata = NULL
  )
)

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

#' Creates a LCMSPlotDataContainer from a features matrix
#'
#' @param raw_data A list of mzR objects
#' @param features_matrix The features matrix
#' @param peaks_matrix The peaks matrix
#' @param samples The samples metadata
#' @param ppm The ppm error for the chromatograms
#' @param rt_tol The RT tolerance for the chromatograms
#' @returns A LCMSPlotDataContainer object
#' @export
setGeneric(
  "create_chromatograms_from_features",
  function(raw_data, features_matrix, peaks_matrix, samples, ppm, rt_tol) standardGeneric("create_chromatograms_from_features")
)

#' @rdname create_chromatograms_from_features
setMethod(
  f = "create_chromatograms_from_features",
  signature = c("list", "data.frame", "data.frame", "data.frame", "numeric", "numeric"),
  definition = function(raw_data, features_matrix, peaks_matrix, samples, ppm, rt_tol) {
    chromatograms <- data.frame()
    mass_traces <- data.frame()
    metadata <- data.frame()

    for (i in seq_len(nrow(features_matrix))) {
      feature <- features_matrix[i,]
      rtr <- c(feature$rt - rt_tol, feature$rt + rt_tol)
      peak_indices <- as.numeric(unlist(strsplit(feature %>% pull(peakidx), ',')))
      peaks <- peaks_matrix %>%
        filter(
          row_number() %in% peak_indices,
          sample %in% samples$sample_index
        )

      for (j in seq_len(nrow(peaks))) {
        peak <- peaks[j,]
        mzdev <- peak$mz * (ppm / 1000000)
        mzr <- c(peak$mz - mzdev, peak$mz + mzdev)
        sample <- samples %>% filter(sample_index == peak$sample)
        raw_obj <- raw_data[[sample$sample_path]]

        data <- .create_chromatogram(raw_obj, mz_range = mzr, rt_range = rtr)
        # raw_eic <- rawEIC(xraw_obj, mzrange = mzr, rtrange = rtr)

        # eic <- data.frame(
        #   intensity = data$chromatograms$intensity,
        #   rt = data$chromatograms$rt, #xraw_obj@scantime[raw_eic$scan],
        #   peak_rt_min = peak$rtmin,
        #   peak_rt_max = peak$rtmax,
        #   sample_id = sample$sample_id,
        #   feature_id = feature$name
        #   # TODO: add sample metadata
        # )

        metadata <- rbind(metadata, data.frame(
          sample_id = sample$sample_id,
          feature_id = feature$name,
          peak_rt_min = peak$rtmin,
          peak_rt_max = peak$rtmax
        ))

        chromatograms <- rbind(chromatograms, data.frame(
          rt = data$chromatograms$rt, #xraw_obj@scantime[raw_eic$scan],
          intensity = data$chromatograms$intensity,
          metadata_index = nrow(metadata)
        ))

        mass_traces <- rbind(mass_traces, data.frame(
          rt = data$mass_traces$rt,
          mz = data$mass_traces$mz,
          metadata_index = nrow(metadata)
        ))
      }
    }

    peakly_chromatogram_container(chromatograms, mass_traces, metadata)
  }
)

#' Creates a LCMSPlotDataContainer from raw sample files
#'
#' @param raw_data A list of mzR objects
#' @param features A matrix of entries with mz and rt values
#' @param samples The samples metadata
#' @param ppm The ppm error for the chromatograms
#' @param rt_tol The RT tolerance for the chromatograms
#' @returns A LCMSPlotDataContainer object
#' @export
setGeneric(
  "create_chromatograms_from_raw",
  function(raw_data, features, samples, ppm, rt_tol) standardGeneric("create_chromatograms_from_raw")
)

#' @rdname create_chromatograms_from_raw
setMethod(
  f = "create_chromatograms_from_raw",
  signature = c("list", "matrix", "data.frame", "numeric", "numeric"),
  definition = function(raw_data, features, samples, ppm, rt_tol) {
    chromatograms <- data.frame()
    mass_traces <- data.frame()
    metadata <- data.frame()
    # chromatograms <- data.frame(rt = numeric(), intensity = numeric())

    for (i in 1:nrow(features)) {
      feature <- features[i,]
      feature_id <- paste0('M', round(feature['mz']), 'T', round(feature['rt']))
      mzdev <- feature['mz'] * (ppm / 1000000)
      mzr <- c(feature['mz'] - mzdev, feature['mz'] + mzdev)
      rtr <- c(feature['rt'] - rt_tol, feature['rt'] + rt_tol)

      for (j in 1:nrow(samples)) {
        sample <- samples[j,]
        raw_obj <- raw_data[[sample$sample_path]]

        data <- .create_chromatogram(raw_obj, mz_range = mzr, rt_range = rtr)

        metadata <- rbind(metadata, data.frame(
          sample_id = sample$sample_id,
          feature_id = feature_id
        ))

        chromatograms <- rbind(chromatograms, data.frame(
          rt = data$chromatograms$rt, #xraw_obj@scantime[raw_eic$scan],
          intensity = data$chromatograms$intensity,
          metadata_index = nrow(metadata)
        ))

        mass_traces <- rbind(mass_traces, data.frame(
          rt = data$mass_traces$rt,
          mz = data$mass_traces$mz,
          metadata_index = nrow(metadata)
        ))

        # chromatograms <- rbind(
        #   chromatograms,
        #   .create_chromatogram(raw_obj, mz_range = mzr, rt_range = rtr) %>%
        #     mutate(sample_id = sample$sample_id,
        #            feature_id = feature_id)
        # )
      }
    }

    peakly_chromatogram_container(chromatograms, mass_traces, metadata)
  }
)

# setAs(
#   from = "LCMSPlotDataContainer",
#   to = "data.frame",
#   def = function(from) {
#     data <- data.frame(
#       intensity = from@intensities,
#       rt = from@rts,
#       sample_id = from@sample_ids,
#       stringsAsFactors = FALSE
#     )
#
#     if (!is.null(from@metadata) & nrow(from@metadata) > 0) {
#       data <- cbind(data, from@metadata)
#     }
#
#     return(data)
#   }
# )
