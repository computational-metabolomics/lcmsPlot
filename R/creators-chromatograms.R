#' Create a lcmsPlotDataContainer from feature IDs.
#'
#' @param obj A lcmsPlotDataContainer object.
#' @param options The plot object's options.
#' @returns A lcmsPlotDataContainer object.
#' @export
#' @examples
#' data_obj <- get_XCMSnExp_object_example(indices = c(1, 7), should_group_peaks = TRUE)
#' data_container <- create_data_container_from_obj(data_obj, sample_id_column = "sample_name", metadata = NULL)
#'
#' opts <- default_options()
#' opts$chromatograms$features <- c('M205T2785', 'M207T2713')
#' opts$chromatograms$sample_ids <- c('ko15', 'wt15')
#' opts$chromatograms$ppm <- 5
#' opts$chromatograms$rt_tol <- 10
#'
#' data_container <- create_chromatograms_from_feature_ids(data_container, opts)
setGeneric(
  "create_chromatograms_from_feature_ids",
  function(obj, options) standardGeneric("create_chromatograms_from_feature_ids")
)

#' @rdname create_chromatograms_from_feature_ids
setMethod(
  f = "create_chromatograms_from_feature_ids",
  signature = c("lcmsPlotDataContainer", "list"),
  definition = function(obj, options) {
    if (!is_xcms_object(obj)) {
      stop("To use feature IDs from the grouped peaks you need to provide an xcms object.")
    }

    metadata <- obj@metadata %>%
      filter(.data$sample_id %in% options$chromatograms$sample_ids)
    raw_data <- io_get_raw_data(metadata$sample_path)
    all_detected_peaks <- get_detected_peaks(obj@data_obj)
    grouped_peaks <- get_grouped_peaks(obj@data_obj) %>%
      filter(.data$name %in% options$chromatograms$features)
    detected_peaks <- data.frame()

    chromatograms <- data.frame()
    mass_traces <- data.frame()
    additional_metadata <- data.frame()

    for (i in seq_len(nrow(grouped_peaks))) {
      feature <- grouped_peaks[i,]
      rtr <- c(
        feature$rt - options$chromatograms$rt_tol,
        feature$rt + options$chromatograms$rt_tol
      )
      peak_indices <- feature %>%
        pull(.data$peakidx) %>%
        strsplit(',') %>%
        unlist() %>%
        as.numeric()
      peaks <- all_detected_peaks %>%
        filter(
          row_number() %in% peak_indices,
          .data$sample_index %in% metadata$sample_index
        ) %>%
        left_join(metadata, by = "sample_index")

      detected_peaks <- rbind(detected_peaks, peaks)

      for (j in seq_len(nrow(peaks))) {
        peak <- peaks[j,]
        mzr <- get_mz_range(peak$mz, options$chromatograms$ppm)
        sample_metadata <- metadata %>%
          filter(.data$sample_index == peak$sample_index)
        raw_obj <- raw_data[[sample_metadata$sample_path]]

        data <- create_chromatogram(
          raw_obj,
          mz_range = mzr,
          rt_range = rtr,
          fill_gaps = options$chromatograms$fill_gaps
        )

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

    io_close_raw_data(raw_data)

    obj@chromatograms <- chromatograms
    obj@mass_traces <- mass_traces
    obj@additional_metadata <- additional_metadata
    obj@detected_peaks <- detected_peaks

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
#' @examples
#' data_obj <- get_XCMSnExp_object_example(indices = c(1, 7))
#' data_container <- create_data_container_from_obj(data_obj, sample_id_column = "sample_name", metadata = NULL)
#'
#' opts <- default_options()
#' opts$chromatograms$features <- rbind(c(mzmin = 334.9, mzmax = 335.1, rtmin = 2700, rtmax = 2750))
#' opts$chromatograms$sample_ids <- c('ko15', 'wt15')
#'
#' data_container <- create_chromatograms_from_features(data_container, opts)
setGeneric(
  "create_chromatograms_from_features",
  function(obj, options) standardGeneric("create_chromatograms_from_features")
)

#' @rdname create_chromatograms_from_features
setMethod(
  f = "create_chromatograms_from_features",
  signature = c("lcmsPlotDataContainer", "list"),
  definition = function(obj, options) {
    metadata <- obj@metadata %>%
      filter(.data$sample_id %in% options$chromatograms$sample_ids)

    raw_data <- io_get_raw_data(metadata$sample_path)
    all_detected_peaks <- get_detected_peaks(obj@data_obj)

    process_sample <- function(i) {
      sample_metadata <- metadata[i, ]
      raw_obj <- raw_data[[sample_metadata$sample_path]]
      hdr <- mzR::header(raw_obj)
      full_rt_range <- range(hdr$retentionTime)

      chromatograms_list <- list()
      mass_traces_list <- list()
      additional_metadata_list <- list()
      detected_peaks_list <- list()

      features <- get_features(
        options,
        sample_metadata,
        full_rt_range = full_rt_range
      )

      for (j in seq_len(length(features))) {
        feature_data <- features[[j]]

        data <- create_chromatogram(
          raw_obj,
          mz_range = feature_data$mzr,
          rt_range = feature_data$rtr,
          fill_gaps = options$chromatograms$fill_gaps
        )

        if (!is.null(all_detected_peaks)) {
          peaks <- all_detected_peaks %>%
            filter(
              .data$sample_index %in% sample_metadata$sample_index,
              .data$mz >= feature_data$mzr[1],
              .data$mz <= feature_data$mzr[2],
              .data$rt >= feature_data$rtr[1],
              .data$rt <= feature_data$rtr[2]
            ) %>%
            left_join(metadata, by = "sample_index")
        } else {
          peaks <- data.frame()
        }

        detected_peaks_list[[j]] <- peaks

        additional_metadata_index <- (i - 1) * nrow(options$chromatograms$features) + j

        additional_metadata_list[[j]] <- data.frame(
          metadata_index = sample_metadata$sample_index,
          feature_id = feature_data$feature_id
        )

        chromatograms_list[[j]] <- data.frame(
          rt = data$chromatograms$rt,
          intensity = data$chromatograms$intensity,
          metadata_index = sample_metadata$sample_index,
          additional_metadata_index = additional_metadata_index
        )

        mass_traces_list[[j]] <- data.frame(
          rt = data$mass_traces$rt,
          mz = data$mass_traces$mz,
          metadata_index = sample_metadata$sample_index,
          additional_metadata_index = additional_metadata_index
        )
      }

      list(
        chromatograms = do.call(rbind, chromatograms_list),
        mass_traces = do.call(rbind, mass_traces_list),
        additional_metadata = do.call(rbind, additional_metadata_list),
        detected_peaks = do.call(rbind, detected_peaks_list)
      )
    }

    if (!is.null(options$parallel_param)) {
      results <- BiocParallel::bplapply(
        seq_len(nrow(metadata)),
        process_sample,
        BPPARAM = options$parallel_param
      )
    } else {
      results <- lapply(seq_len(nrow(metadata)), process_sample)
    }

    io_close_raw_data(raw_data)

    obj@chromatograms <- do.call(rbind, lapply(results, `[[`, "chromatograms"))
    obj@mass_traces <- do.call(rbind, lapply(results, `[[`, "mass_traces"))
    obj@additional_metadata <- do.call(rbind, lapply(results, `[[`, "additional_metadata"))
    obj@detected_peaks <- do.call(rbind, lapply(results, `[[`, "detected_peaks"))

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
#' @examples
#' data_obj <- get_XCMSnExp_object_example(indices = c(1, 7))
#' data_container <- create_data_container_from_obj(data_obj, sample_id_column = "sample_name", metadata = NULL)
#'
#' opts <- default_options()
#' opts$chromatograms$aggregation_fun <- 'max'
#' opts$chromatograms$sample_ids <- c('ko15', 'wt15')
#'
#' data_container <- create_full_rt_chromatograms(data_container, opts)
setGeneric(
  "create_full_rt_chromatograms",
  function(obj, options) standardGeneric("create_full_rt_chromatograms")
)

#' @rdname create_full_rt_chromatograms
setMethod(
  f = "create_full_rt_chromatograms",
  signature = c("lcmsPlotDataContainer", "list"),
  definition = function(obj, options) {
    metadata <- obj@metadata %>%
      filter(.data$sample_id %in% options$chromatograms$sample_ids)
    raw_data <- io_get_raw_data(metadata$sample_path)

    process_sample <- function(i) {
      sample_metadata <- metadata[i, ]
      raw_obj <- raw_data[[sample_metadata$sample_path]]

      if (options$chromatograms$rt_adjusted) {
        if (!is_xcms_object(obj)) {
          stop("The data object should be XCMSnExp or MsExperiment to plot the RT adjusted chromatograms")
        }

        scan_indices <- which(xcms::fromFile(obj@data_obj) == sample_metadata$sample_index)
        rt_adjusted <- xcms::adjustedRtime(obj@data_obj)[scan_indices]
      } else {
        rt_adjusted <- NULL
      }

      data <- create_bpc_tic(
        raw_obj,
        options$chromatograms$aggregation_fun,
        rt_adjusted
      )

      data.frame(
        rt = data$chromatograms$rt,
        intensity = data$chromatograms$intensity,
        metadata_index = sample_metadata$sample_index,
        additional_metadata_index = sample_metadata$sample_index
      )
    }

    if (!is.null(options$parallel_param)) {
      chromatograms_list <- BiocParallel::bplapply(
        seq_len(nrow(metadata)),
        process_sample,
        BPPARAM = options$parallel_param)
    } else {
      chromatograms_list <- lapply(seq_len(nrow(metadata)), process_sample)
    }

    chromatograms <- do.call(rbind, chromatograms_list)

    io_close_raw_data(raw_data)

    obj@chromatograms <- chromatograms

    validObject(obj)

    return(obj)
  }
)
