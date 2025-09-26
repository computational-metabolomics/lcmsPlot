#' Creates/updates an lcmsPlotDataContainer from spectra
#'
#' @param obj A lcmsPlotDataContainer object
#' @param options The lcmsPlot options
#' @returns A lcmsPlotDataContainer object
#' @export
#' @examples
#' data_obj <- get_XCMSnExp_object_example(indices = 1)
#' data_container <- create_data_container_from_obj(data_obj, sample_id_column = "sample_name", metadata = NULL)
#'
#' opts <- default_options()
#' opts$spectra$show <- TRUE
#' opts$spectra$sample_ids <- 'ko15'
#' opts$spectra$scan_index <- 456
#'
#' data_container <- create_spectra(data_container, opts)
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

    metadata <- obj@metadata %>%
      filter(.data$sample_id %in% options$spectra$sample_ids)

    grouped_peaks <- get_grouped_peaks(obj@data_obj) %>%
      { if (!is.null(.) && nrow(.) > 0) filter(., .data$name %in% options$chromatograms$features) else NULL } %>%
      { if (!is.null(.) && nrow(.) == 0) NULL else . }

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

    # Process each sample
    for (i in seq_len(nrow(metadata))) {
      sample_metadata <- metadata[i,]
      raw_obj <- mzR::openMSfile(sample_metadata$sample_path)

      if (is_standalone) {
        spectra <- create_spectra_for_sample(
          raw_obj,
          obj@detected_peaks,
          sample_metadata,
          options)

        spectra <- spectra %>%
          mutate(
            metadata_index = sample_metadata$sample_index,
            additional_metadata_index = NA
          )
        all_spectra <- rbind(all_spectra, spectra)
      } else {
        features <- get_features(options, sample_metadata, grouped_peaks = grouped_peaks)

        for (j in seq_len(length(features))) {
          feature_data <- features[[j]]

          spectra <- create_spectra_for_sample(
            raw_obj,
            obj@detected_peaks,
            sample_metadata,
            options,
            rt_range = feature_data$rtr)

          spectra <- spectra %>%
            mutate(
              metadata_index = sample_metadata$sample_index,
              additional_metadata_index = additional_metadata_index
            )
          all_spectra <- rbind(all_spectra, spectra)

          additional_metadata_index <- additional_metadata_index + 1
        }
      }

      mzR::close(raw_obj)
    }

    all_spectra <- all_spectra %>%
      mutate(reference = FALSE)

    if (!is.null(spectral_library)) {
      all_spectra_as_list <- all_spectra %>%
        group_by(.data$metadata_index, .data$rt) %>%
        group_split()

      query_spectra <- S4Vectors::DataFrame(
        metadata_index = unique(all_spectra$metadata_index),
        rt = unique(all_spectra$rt)
      )
      query_spectra$mz <- lapply(all_spectra_as_list, function (x) x$mz)
      query_spectra$intensity <- lapply(all_spectra_as_list, function (x) x$intensity)
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

      i <- 1
      for (target_index in target_indices) {
        hit_mzs <- Spectra::mz(spectral_library)[[target_index]]
        hit_intensities <- Spectra::intensity(spectral_library)[[target_index]]

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
      group_by(.data$metadata_index, .data$rt, .data$reference) %>%
      mutate(intensity = 100 * .data$intensity / max(.data$intensity)) %>%
      ungroup()

    obj@spectra <- all_spectra

    validObject(obj)

    return(obj)
  }
)
