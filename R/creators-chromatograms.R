#' Create an instance of class `lcmsPlotDataContainer` from feature IDs
#'
#' The input features are specified in the `options` list
#' under `chromatograms$features`.
#'
#' @param obj An instance of class `lcmsPlotDataContainer`.
#' @param options A `list` representing the plot object's options.
#' @return An instance of class `lcmsPlotDataContainer` with chromatograms
#' of the specified feature IDs.
#' @keywords internal
setGeneric(
    "create_chromatograms_from_feature_ids",
    function(obj, options)
        standardGeneric("create_chromatograms_from_feature_ids")
)

#' @rdname create_chromatograms_from_feature_ids
setMethod(
    f = "create_chromatograms_from_feature_ids",
    signature = c("lcmsPlotDataContainer", "list"),
    definition = function(obj, options) {
        if (!is_xcms_data(obj@data_obj)) {
            stop("To use feature IDs from the grouped peaks you need to provide an xcms object.")
        }

        metadata <- obj@metadata |>
            filter(.data$sample_id %in% options$chromatograms$sample_ids)
        raw_data <- io_get_raw_data(metadata$sample_path)
        adjusted_rts <- get_adjusted_rts(obj@data_obj)
        all_detected_peaks <- get_detected_peaks(obj@data_obj)
        grouped_peaks <- get_grouped_peaks(obj@data_obj) |>
            filter(.data$name %in% options$chromatograms$features)
        detected_peaks <- data.frame()

        chromatograms <- data.frame()
        mass_traces <- data.frame(
            rt = numeric(),
            mz = numeric(),
            metadata_index = numeric(),
            additional_metadata_index = numeric()
        )
        additional_metadata <- data.frame()

        for (i in seq_len(nrow(grouped_peaks))) {
            feature <- grouped_peaks[i,]
            rtr <- c(
                feature$rt - options$chromatograms$rt_tol,
                feature$rt + options$chromatograms$rt_tol
            )
            peak_indices <- feature |>
                pull(.data$peakidx) |>
                strsplit(',') |>
                unlist() |>
                as.numeric()
            peaks <- all_detected_peaks |>
                filter(
                    row_number() %in% peak_indices,
                    .data$sample_index %in% metadata$sample_index
                ) |>
                left_join(metadata, by = "sample_index")

            detected_peaks <- rbind(detected_peaks, peaks)

            for (j in seq_len(nrow(peaks))) {
                peak <- peaks[j,]
                mzr <- get_mz_range(peak$mz, options$chromatograms$ppm)
                sample_metadata <- metadata |>
                    filter(.data$sample_index == peak$sample_index)
                raw_obj <- raw_data[[sample_metadata$sample_path]]

                if (!is.null(adjusted_rts)) {
                    sample_adjusted_rt <- adjusted_rts |>
                        filter(.data$file_index == sample_metadata$sample_index)
                } else {
                    sample_adjusted_rt <- NULL
                }

                data <- create_chromatogram(
                    raw_obj,
                    mz_range = mzr,
                    rt_range = rtr,
                    fill_gaps = options$chromatograms$fill_gaps,
                    adjusted_rt = sample_adjusted_rt
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

                if (nrow(data$mass_traces) > 0) {
                    mass_traces <- rbind(mass_traces, data.frame(
                        rt = data$mass_traces$rt,
                        mz = data$mass_traces$mz,
                        metadata_index = sample_metadata$sample_index,
                        additional_metadata_index = nrow(additional_metadata)
                    ))
                }
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

#' Creates an instance of class `lcmsPlotDataContainer` from a features `matrix`
#' or `data.frame`
#'
#' The input features are specified in the `options` list
#' under `chromatograms$features`.
#'
#' @param obj An instance of class `lcmsPlotDataContainer`.
#' @param options A `list` representing the plot object's options.
#' @return An instance of class `lcmsPlotDataContainer` with chromatograms
#' of the specified features.
#' @keywords internal
setGeneric(
    "create_chromatograms_from_features",
    function(obj, options) standardGeneric("create_chromatograms_from_features")
)

#' @rdname create_chromatograms_from_features
setMethod(
    f = "create_chromatograms_from_features",
    signature = c("lcmsPlotDataContainer", "list"),
    definition = function(obj, options) {
        metadata <- obj@metadata |>
            filter(.data$sample_id %in% options$chromatograms$sample_ids)

        raw_data <- io_get_raw_data(metadata$sample_path)
        all_detected_peaks <- get_detected_peaks(obj@data_obj)
        adjusted_rts <- get_adjusted_rts(obj@data_obj)

        process_sample <- function(i) {
            sample_metadata <- metadata[i, ]
            raw_obj <- raw_data[[sample_metadata$sample_path]]

            if (!is.null(adjusted_rts)) {
                sample_adjusted_rt <- adjusted_rts |>
                    filter(.data$file_index == sample_metadata$sample_index)
            } else {
                sample_adjusted_rt <- NULL
            }

            hdr <- ms_header(raw_obj)
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
                    fill_gaps = options$chromatograms$fill_gaps,
                    adjusted_rt = sample_adjusted_rt
                )

                if (!is.null(all_detected_peaks)) {
                    peaks <- all_detected_peaks |>
                        filter(
                            .data$sample_index %in% sample_metadata$sample_index,
                            .data$mz >= feature_data$mzr[1],
                            .data$mz <= feature_data$mzr[2],
                            .data$rt >= feature_data$rtr[1],
                            .data$rt <= feature_data$rtr[2]
                        ) |>
                        left_join(metadata, by = "sample_index")
                } else {
                    peaks <- data.frame()
                }

                detected_peaks_list[[j]] <- peaks

                n_features <- nrow(options$chromatograms$features)
                additional_metadata_index <- (i - 1) * n_features + j

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

                if (nrow(data$mass_traces) > 0) {
                    mass_traces_list[[j]] <- data.frame(
                        rt = data$mass_traces$rt,
                        mz = data$mass_traces$mz,
                        metadata_index = sample_metadata$sample_index,
                        additional_metadata_index = additional_metadata_index
                    )
                }
            }

            if (length(mass_traces_list) == 0) {
                mass_traces <- data.frame(
                    rt = numeric(),
                    mz = numeric(),
                    metadata_index = numeric(),
                    additional_metadata_index = numeric()
                )
            } else {
                mass_traces <- do.call(rbind, mass_traces_list)
            }

            list(
                chromatograms = do.call(rbind, chromatograms_list),
                mass_traces = mass_traces,
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

        obj@chromatograms <- do.call(
            rbind,
            lapply(results, `[[`, "chromatograms"))
        obj@mass_traces <- do.call(
            rbind,
            lapply(results, `[[`, "mass_traces"))
        obj@additional_metadata <- do.call(
            rbind,
            lapply(results, `[[`, "additional_metadata"))
        obj@detected_peaks <- do.call(
            rbind,
            lapply(results, `[[`, "detected_peaks"))

        validObject(obj)

        return(obj)
    }
)

#' Creates an instance of class `lcmsPlotDataContainer` from
#' base peak chromatograms (BPC) or total ion chromatograms (TIC)
#'
#' The type of chromatograms is defined in
#' `options$chromatograms$aggregation_fun`.
#'
#' @param obj An instance of class `lcmsPlotDataContainer`.
#' @param options A `list` representing the plot object's options.
#' @return An instance of class `lcmsPlotDataContainer` with BPC or TIC.
#' @keywords internal
setGeneric(
    "create_full_rt_chromatograms",
    function(obj, options) standardGeneric("create_full_rt_chromatograms")
)

#' @rdname create_full_rt_chromatograms
setMethod(
    f = "create_full_rt_chromatograms",
    signature = c("lcmsPlotDataContainer", "list"),
    definition = function(obj, options) {
        metadata <- obj@metadata |>
            filter(.data$sample_id %in% options$chromatograms$sample_ids)
        raw_data <- io_get_raw_data(metadata$sample_path)

        process_sample <- function(i) {
            sample_metadata <- metadata[i, ]
            raw_obj <- raw_data[[sample_metadata$sample_path]]

            if (options$chromatograms$rt_type %in% c("corrected", "both")) {
                if (!is_xcms_data(obj@data_obj)) {
                    stop("The data object should be XCMSnExp or MsExperiment to plot the RT adjusted chromatograms")
                }

                scan_indices <- which(
                    xcms::fromFile(obj@data_obj) == sample_metadata$sample_index
                )
                rt_adjusted <- xcms::adjustedRtime(obj@data_obj)[scan_indices]
            } else {
                rt_adjusted <- NULL
            }

            data <- create_bpc_tic(
                raw_obj,
                options$chromatograms$aggregation_fun
            )
            chroms <- data$chromatograms

            if (options$chromatograms$rt_type == "both") {
                additional_metadata_index <- c(
                    rep(1, nrow(chroms)),
                    rep(2, nrow(chroms))
                )
                chroms_rt_adjusted <- chroms
                chroms_rt_adjusted$rt <- rt_adjusted
                chroms <- rbind(chroms, chroms_rt_adjusted)
            } else if (options$chromatograms$rt_type == "corrected") {
                chroms$rt <- rt_adjusted
                additional_metadata_index <- rep(2, nrow(chroms))
            } else {
                additional_metadata_index <- rep(1, nrow(chroms))
            }

            chromatograms <- data.frame(
                rt = chroms$rt,
                intensity = chroms$intensity,
                metadata_index = sample_metadata$sample_index,
                additional_metadata_index = additional_metadata_index
            )
        }

        if (!is.null(options$parallel_param)) {
            chromatograms_list <- BiocParallel::bplapply(
                seq_len(nrow(metadata)),
                process_sample,
                BPPARAM = options$parallel_param)
        } else {
            chromatograms_list <- lapply(
                seq_len(nrow(metadata)),
                process_sample)
        }

        chromatograms <- do.call(rbind, chromatograms_list)

        io_close_raw_data(raw_data)

        obj@chromatograms <- chromatograms
        obj@additional_metadata <- data.frame(
            rt_adjusted = c("RT uncorrected", "RT corrected")
        )

        validObject(obj)

        return(obj)
    }
)

#' Creates an instance of class `lcmsPlotDataContainer` from a
#' Compound Discoverer results object.
#'
#' @param obj An instance of class `lcmsPlotDataContainer`.
#' @param options A `list` representing the plot object's options.
#' @return An instance of class `lcmsPlotDataContainer` with chromatograms
#' of the specified Compound Discoverer compounds.
#' @keywords internal
setGeneric(
    "create_chromatograms_from_compound_discoverer",
    function(obj, options)
        standardGeneric("create_chromatograms_from_compound_discoverer")
)

#' @rdname create_chromatograms_from_compound_discoverer
setMethod(
    f = "create_chromatograms_from_compound_discoverer",
    signature = c("lcmsPlotDataContainer", "list"),
    definition = function(obj, options) {
        cd_opts <- options$compound_discoverer

        metadata <- obj@metadata |>
            filter(.data$sample_id %in% options$chromatograms$sample_ids)

        xic_traces_results <- get_xic_traces_from_compounds(
            obj@data_obj,
            cd_opts$compounds_query)

        compounds <- xic_traces_results |>
            group_by(.data$name) |>
            distinct(
                .data$name,
                .data$formula,
                .data$mz,
                .data$adduct
            ) |>
            ungroup() |>
            as.data.frame() |>
            mutate(index = row_number())

        process_sample <- function(i) {
            sample_metadata <- metadata[i, ]

            chromatograms_list <- list()
            additional_metadata_list <- list()
            detected_peaks_list <- list()

            for (j in seq_len(nrow(compounds))) {
                compound_data <- compounds[j,]
                cols <- names(compound_data)
                cols <- cols[cols != "index"]

                xic_entry <- xic_traces_results |>
                    filter(.data$sample_id == sample_metadata$sample_id) |>
                    filter(
                        across(
                            all_of(cols),
                            ~ . == compound_data[[cur_column()]]
                        )
                    ) |>
                    as.data.frame()

                if (nrow(xic_entry) > 0) {
                    chroms <- parse_trace(xic_entry$trace[[1]]) |>
                        filter(
                            .data$rt >= xic_entry$rtmin - cd_opts$rt_extend,
                            .data$rt <= xic_entry$rtmax + cd_opts$rt_extend)

                    chromatograms_list[[j]] <- data.frame(
                        rt = chroms$rt,
                        intensity = chroms$intensity,
                        metadata_index = sample_metadata$sample_index,
                        additional_metadata_index = compound_data$index
                    )

                    additional_metadata_list[[j]] <- data.frame(
                        metadata_index = sample_metadata$sample_index,
                        name = compound_data$name,
                        formula = compound_data$formula,
                        mz = compound_data$mz,
                        adduct = compound_data$adduct
                    )

                    detected_peaks_list[[j]] <- xic_entry |>
                        select(
                            .data$name,
                            .data$sample_id,
                            .data$mz,
                            .data$rt,
                            .data$rtmin,
                            .data$rtmax,
                            .data$into,
                            .data$maxo
                        ) |>
                        left_join(metadata, by = "sample_id")
                }
            }

            list(
                chromatograms = do.call(rbind, chromatograms_list),
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

        obj@chromatograms <- do.call(
            rbind,
            lapply(results, `[[`, "chromatograms"))
        obj@additional_metadata <- do.call(
            rbind,
            lapply(results, `[[`, "additional_metadata"))
        obj@detected_peaks <- do.call(
            rbind,
            lapply(results, `[[`, "detected_peaks"))

        validObject(obj)

        return(obj)
    }
)
