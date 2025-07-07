.validators <- list(
  additional_metadata = function(df) {
    nrow(df) == 0 || "metadata_index" %in% colnames(df)
  },
  chromatograms = function(df) {
    nrow(df) == 0 || identical(colnames(df), c("rt", "intensity", "metadata_index", "additional_metadata_index"))
  },
  mass_traces = function(df) {
    nrow(df) == 0 || identical(colnames(df), c("rt", "mz", "metadata_index", "additional_metadata_index"))
  },
  spectra = function(df) {
    nrow(df) == 0 || identical(colnames(df), c("mz", "intensity", "rt", "metadata_index", "additional_metadata_index", "reference"))
  },
  intensity_maps = function(df) {
    nrow(df) == 0 || identical(colnames(df), c("rt", "mz", "intensity", "metadata_index", "additional_metadata_index"))
  },
  detected_peaks = function(df) {
    nrow(df) == 0 || all(c("mz", "rt", "rtmin", "rtmax", "sample_index") %in% colnames(df))
  }
)

#' Create an lcmsPlotDataContainer object from a data object (e.g., XCMSnExp)
#'
#' @param data_obj The data object (e.g., XCMSnExp)
#' @param sample_id_column The sample ID column
#' @param metadata The sample metadata in case it's not provided in the data object
create_data_container_from_obj <- function(data_obj, sample_id_column, metadata) {
  new("lcmsPlotDataContainer",
      data_obj = data_obj,
      metadata = get_metadata(data_obj, sample_id_column, metadata),
      chromatograms = data.frame(),
      mass_traces = data.frame(),
      spectra = data.frame(),
      intensity_maps = data.frame(),
      additional_metadata = data.frame(),
      detected_peaks = data.frame())
}

#' lcmsPlotDataContainer class
#'
#' @slot data_obj The data object. One of:
#' - XCMSnExp
#' - MsExperiment
#' - character: vector of sample paths
#' @slot metadata The sample metadata
#' @slot chromatograms The chromatograms
#' @slot mass_traces The mass traces
#' @slot spectra The spectra
#' @slot intensity_maps The 2D intensity maps
#' @slot additional_metadata Additional information attached to datasets
#' @slot detected_peaks The detected peaks
#' @export
setClass(
  "lcmsPlotDataContainer",
  slots = list(
    data_obj = "ANY",
    metadata = "data.frame",
    chromatograms = "data.frame",
    mass_traces = "data.frame",
    spectra = "data.frame",
    intensity_maps = "data.frame",
    additional_metadata = "data.frame",
    detected_peaks = "data.frame"
  ),
  prototype = list(
    data_obj = NULL,
    metadata = NULL,
    chromatograms = NULL,
    mass_traces = NULL,
    spectra = NULL,
    intensity_maps = NULL,
    additional_metadata = NULL,
    detected_peaks = NULL
  )
)

setValidity("lcmsPlotDataContainer", function(object) {
  ret <- TRUE

  if (!inherits(object@data_obj, c("XCMSnExp", "MsExperiment", "character"))) {
    ret <- "@data_obj must be either XCMSnExp, MsExperiment, or character."
  } else {
    for (validator_name in names(.validators)) {
      df <- slot(object, validator_name)

      if (!.validators[[validator_name]](df)) {
        ret <- paste0(validator_name, " did not pass validation.")
        break
      }
    }
  }

  ret
})

#' Creates a lcmsPlotDataContainer from feature IDs
#'
#' @param obj A lcmsPlotDataContainer object
#' @param options The options
#' @returns A lcmsPlotDataContainer object
#' @export
setGeneric(
  "create_chromatograms_from_feature_ids",
  function(obj, options) standardGeneric("create_chromatograms_from_feature_ids")
)

#' @rdname create_chromatograms_from_feature_ids
setMethod(
  f = "create_chromatograms_from_feature_ids",
  signature = c("lcmsPlotDataContainer", "list"),
  definition = function(obj, options) {
    metadata <- obj@metadata %>% filter(sample_id %in% options$chromatograms$sample_ids)
    raw_data <- io_utils$get_raw_data(metadata$sample_path)
    all_detected_peaks <- get_detected_peaks(obj@data_obj)
    grouped_peaks <- get_grouped_peaks(obj@data_obj) %>% filter(name %in% options$chromatograms$features)
    detected_peaks <- data.frame()

    chromatograms <- data.frame()
    mass_traces <- data.frame()
    additional_metadata <- data.frame()

    for (i in seq_len(nrow(grouped_peaks))) {
      feature <- grouped_peaks[i,]
      rtr <- c(feature$rt - options$chromatograms$rt_tol, feature$rt + options$chromatograms$rt_tol)
      peak_indices <- as.numeric(unlist(strsplit(feature %>% pull(peakidx), ',')))
      peaks <- all_detected_peaks %>%
        filter(
          row_number() %in% peak_indices,
          sample_index %in% metadata$sample_index
        ) %>%
        left_join(metadata, by = "sample_index")

      detected_peaks <- rbind(detected_peaks, peaks)

      for (j in seq_len(nrow(peaks))) {
        peak <- peaks[j,]
        mzr <- get_mz_range(peak$mz, options$chromatograms$ppm)
        sample_metadata <- metadata %>% filter(sample_index == peak$sample_index)
        raw_obj <- raw_data[[sample_metadata$sample_path]]

        data <- create_chromatogram(raw_obj, mz_range = mzr, rt_range = rtr)

        additional_metadata <- rbind(additional_metadata, data.frame(
          metadata_index = sample_metadata$sample_index,
          feature_id = feature$name
        ))

        chromatograms <- rbind(chromatograms, data.frame(
          rt = data$chromatograms$rt,
          intensity = data$chromatograms$intensity,
          metadata_index = sample_metadata$sample_index,
          additional_metadata_index = nrow(additional_metadata)
        ))

        mass_traces <- rbind(mass_traces, data.frame(
          rt = data$mass_traces$rt,
          mz = data$mass_traces$mz,
          metadata_index = sample_metadata$sample_index,
          additional_metadata_index = nrow(additional_metadata)
        ))
      }
    }

    io_utils$close_raw_data(raw_data)

    obj@chromatograms <- chromatograms
    obj@mass_traces <- mass_traces
    obj@additional_metadata <- additional_metadata
    obj@detected_peaks = detected_peaks

    validObject(obj)

    return(obj)
  }
)

#' Creates a lcmsPlotDataContainer from a features matrix
#'
#' @param obj A lcmsPlotDataContainer object
#' @param options The options
#' @returns A lcmsPlotDataContainer object
#' @export
setGeneric(
  "create_chromatograms_from_features",
  function(obj, options) standardGeneric("create_chromatograms_from_features")
)

#' @rdname create_chromatograms_from_features
setMethod(
  f = "create_chromatograms_from_features",
  signature = c("lcmsPlotDataContainer", "list"),
  definition = function(obj, options) {
    metadata <- obj@metadata %>% filter(sample_id %in% options$chromatograms$sample_ids)
    raw_data <- io_utils$get_raw_data(metadata$sample_path)
    all_detected_peaks <- get_detected_peaks(obj@data_obj)

    chromatograms <- data.frame()
    mass_traces <- data.frame()
    additional_metadata <- data.frame()
    detected_peaks <- data.frame()

    for (i in 1:nrow(options$chromatograms$features)) {
      feature <- options$chromatograms$features[i,]
      feature_id <- paste0('M', round(feature['mz']), 'T', round(feature['rt']))
      mzr <- get_mz_range(feature['mz'], options$chromatograms$ppm)
      rtr <- c(feature['rt'] - options$chromatograms$rt_tol, feature['rt'] + options$chromatograms$rt_tol)

      for (j in 1:nrow(metadata)) {
        sample_metadata <- metadata[j,]
        raw_obj <- raw_data[[sample_metadata$sample_path]]

        data <- create_chromatogram(raw_obj, mz_range = mzr, rt_range = rtr)

        if (!is.null(all_detected_peaks)) {
          peaks <- all_detected_peaks %>%
            filter(
              sample_index %in% sample_metadata$sample_index,
              mz >= mzr[1],
              mz <= mzr[2],
              rt >= rtr[1],
              rt <= rtr[2]
            ) %>%
            left_join(metadata, by = "sample_index")
        } else {
          peaks <- data.frame()
        }

        detected_peaks <- rbind(detected_peaks, peaks)

        additional_metadata <- rbind(additional_metadata, data.frame(
          metadata_index = sample_metadata$sample_index,
          feature_id = feature_id
        ))

        chromatograms <- rbind(chromatograms, data.frame(
          rt = data$chromatograms$rt,
          intensity = data$chromatograms$intensity,
          metadata_index = sample_metadata$sample_index,
          additional_metadata_index = nrow(additional_metadata)
        ))

        mass_traces <- rbind(mass_traces, data.frame(
          rt = data$mass_traces$rt,
          mz = data$mass_traces$mz,
          metadata_index = sample_metadata$sample_index,
          additional_metadata_index = nrow(additional_metadata)
        ))
      }
    }

    io_utils$close_raw_data(raw_data)

    obj@chromatograms <- chromatograms
    obj@mass_traces <- mass_traces
    obj@additional_metadata <- additional_metadata
    obj@detected_peaks <- detected_peaks

    validObject(obj)

    return(obj)
  }
)

#' Creates a lcmsPlotDataContainer from BPC or TIC data
#'
#' @param obj A lcmsPlotDataContainer object
#' @param options The options
#' @returns A lcmsPlotDataContainer object
#' @export
setGeneric(
  "create_full_rt_chromatograms",
  function(obj, options) standardGeneric("create_full_rt_chromatograms")
)

#' @rdname create_full_rt_chromatograms
setMethod(
  f = "create_full_rt_chromatograms",
  signature = c("lcmsPlotDataContainer", "list"),
  definition = function(obj, options) {
    metadata <- obj@metadata %>% filter(sample_id %in% options$chromatograms$sample_ids)
    raw_data <- io_utils$get_raw_data(metadata$sample_path)

    chromatograms <- data.frame()

    for (i in seq_len(nrow(metadata))) {
      sample_metadata <- metadata[i,]
      raw_obj <- raw_data[[sample_metadata$sample_path]]

      data <- create_bpc_tic(raw_obj, options$chromatograms$aggregation_fun)

      chromatograms <- rbind(chromatograms, data.frame(
        rt = data$chromatograms$rt,
        intensity = data$chromatograms$intensity,
        metadata_index = i,
        additional_metadata_index = i
      ))
    }

    io_utils$close_raw_data(raw_data)

    obj@chromatograms <- chromatograms

    validObject(obj)

    return(obj)
  }
)

#' Creates/updates an lcmsPlotDataContainer from spectra
#'
#' @param obj A lcmsPlotDataContainer object
#' @param options The lcmsPlot options
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
    is_standalone <- !options$chromatograms$show & options$spectra$show

    if  (is_standalone) {
      metadata <- obj@metadata %>%
        filter(sample_id %in% options$spectra$sample_ids)
    } else {
      metadata <- obj@metadata %>%
        filter(sample_id %in% options$chromatograms$sample_ids)
    }

    if (is_standalone) {
      grouped_peaks <- NULL
    } else {
      grouped_peaks <- get_grouped_peaks(obj@data_obj) %>%
        filter(name %in% options$chromatograms$features)
    }

    all_spectra <- data.frame()

    # Retrieve the spectral library
    spectral_library <- NULL
    if (!is.null(options$spectra$spectral_match_db)) {
      if (length(options$spectra$spectral_match_db) == 1 & all(endsWith(options$spectra$spectral_match_db, ".msp"))) {
        source <- MsBackendMsp::MsBackendMsp()
      } else {
        source <- Spectra::MsBackendMzR()
      }

      spectral_library <- Spectra::Spectra(options$spectra$spectral_match_db, source = source)
    }

    additional_metadata_index <- 1

    for (i in seq_len(nrow(metadata))) {
      sample_metadata <- metadata[i,]
      raw_obj <- mzR::openMSfile(sample_metadata$sample_path)

      if (is_standalone) {
        spectra <- create_spectra_for_sample(
          raw_obj,
          obj@detected_peaks,
          sample_metadata$sample_id,
          options)

        spectra <- spectra %>%
          mutate(
            metadata_index = i,
            additional_metadata_index = NA
          )
        all_spectra <- rbind(all_spectra, spectra)
      } else {
        n_features <- ifelse(
          is.matrix(options$chromatograms$features),
          nrow(options$chromatograms$features),
          length(options$chromatograms$features)
        )

        rt_tol <- options$chromatograms$rt_tol

        for (j in seq_len(n_features)) {
          if (is.matrix(options$chromatograms$features)) {
            feature <- options$chromatograms$features[j,]
            feature_id <- paste0('M', round(feature['mz']), 'T', round(feature['rt']))
            rtr <- c(feature['rt'] - rt_tol, feature['rt'] + rt_tol)
          } else {
            feature <- grouped_peaks[j,]
            feature_id <- feature$name
            rtr <- c(feature$rt - rt_tol, feature$rt + rt_tol)
          }

          spectra <- create_spectra_for_sample(
            raw_obj,
            obj@detected_peaks,
            sample_metadata$sample_id,
            options,
            rt_range = rtr)

          spectra <- spectra %>%
            mutate(
              metadata_index = i,
              additional_metadata_index = additional_metadata_index
            )
          all_spectra <- rbind(all_spectra, spectra)

          additional_metadata_index <- additional_metadata_index + 1
        }
      }

      mzR::close(raw_obj)
    }

    all_spectra = all_spectra %>%
      mutate(reference = FALSE)

    if (!is.null(spectral_library)) {
      all_spectra_as_list <- all_spectra %>%
        group_by(metadata_index, rt) %>%
        group_split()

      query_spectra <- S4Vectors::DataFrame(
        metadata_index = unique(all_spectra$metadata_index),
        rt = unique(all_spectra$rt)
      )
      query_spectra$mz = lapply(all_spectra_as_list, function (x) x$mz)
      query_spectra$intensity = lapply(all_spectra_as_list, function (x) x$intensity)
      query_spectra <- Spectra::Spectra(query_spectra)

      similarities <- Spectra::compareSpectra(query_spectra, spectral_library)

      target_values <- apply(similarities, 1, function(row) {
        sorted_unique <- sort(unique(row), decreasing = TRUE)
        if (length(sorted_unique) >= options$spectra$match_target_index) {
          sorted_unique[options$spectra$match_target_index]
        } else {
          NA
        }
      })

      target_indices <- lapply(1:nrow(similarities), function(i) {
        which(similarities[i, ] == target_values[i])
      }) %>% unlist()

      i = 1
      for (target_index in target_indices) {
        hit_mzs <- mz(spectral_library)[[target_index]]
        hit_intensities <- intensity(spectral_library)[[target_index]]

        all_spectra <- rbind(all_spectra, data.frame(
          mz = hit_mzs,
          intensity = hit_intensities,
          rt = unique(all_spectra_as_list[[i]]$rt),
          metadata_index = unique(all_spectra_as_list[[i]]$metadata_index),
          additional_metadata_index = unique(all_spectra_as_list[[i]]$additional_metadata_index),
          reference = TRUE
        ))

        i <- i + 1
      }
    }

    all_spectra <- all_spectra %>%
      group_by(metadata_index, rt, reference) %>%
      mutate(intensity = 100 * intensity / max(intensity)) %>%
      ungroup()

    obj@spectra <- all_spectra

    validObject(obj)

    return(obj)
  }
)

#' Creates an lcmsPlotDataContainer from an intensity map
#'
#' @param obj A lcmsPlotDataContainer object
#' @param options The lcmsPlot options
#' @returns A lcmsPlotDataContainer object
#' @export
setGeneric(
  "create_intensity_map",
  function(obj, options) standardGeneric("create_intensity_map")
)

#' @rdname create_intensity_map
setMethod(
  f = "create_intensity_map",
  signature = c("lcmsPlotDataContainer", "list"),
  definition = function(obj, options) {
    intensity_maps <- data.frame()

    for (i in seq_len(nrow(obj@metadata))) {
      sample_metadata <- obj@metadata[i,]
      ms <- mzR::openMSfile(sample_metadata$sample_path)
      hdr <- mzR::header(ms)

      rt_range <- options$intensity_maps$rt_range
      mz_range <- options$intensity_maps$mz_range

      idx <- which(hdr$retentionTime >= rt_range[1] & hdr$retentionTime <= rt_range[2])

      scans <- lapply(idx, function(i) {
        pk <- mzR::peaks(ms, i)
        if (nrow(pk) > 0) {
          pk <- pk[pk[,1] >= mz_range[1] & pk[,1] <= mz_range[2], ]
          if (nrow(pk) > 0) {
            data.frame(rt = rep(hdr$retentionTime[i], nrow(pk)),
                       mz = pk[,1],
                       intensity = pk[,2])
          }
        }
      })
      df <- do.call(rbind, scans)

      mzR::close(ms)

      intensity_map <- df %>%
        mutate(rt = round(rt, 1),
               mz = round(mz, 1)) %>%
        group_by(rt, mz) %>%
        summarize(intensity = sum(intensity), .groups = 'drop') %>%
        mutate(metadata_index = i, additional_metadata_index = i)

      intensity_maps <- rbind(intensity_maps, intensity_map)
    }

    obj@intensity_maps <- intensity_maps

    validObject(obj)

    return(obj)
  }
)
