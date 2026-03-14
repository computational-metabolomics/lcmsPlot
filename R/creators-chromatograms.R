
#' Create chromatogram data from a data object
#'
#' Dispatches to the appropriate implementation based on the type of `data_obj`
#' and `features`. Returns a named list with slots to be assigned onto
#' `lcmsPlotDataContainer`.
#'
#' @param data_obj The data object (e.g. `XCMSnExp`, `DBIConnection`, `character`).
#' @param metadata A `data.frame` of sample metadata.
#' @param options A `list` of plot options.
#' @param features `NULL`, a `character` vector of feature IDs, or a
#' `matrix`/`data.frame` of feature ranges.
#' @return A named `list` with elements `chromatograms`, `mass_traces`,
#' `feature_metadata`, and `detected_peaks`.
#' @keywords internal
setGeneric(
    "create_chromatograms",
    function(data_obj, metadata, options, features)
        standardGeneric("create_chromatograms")
)

#' @rdname create_chromatograms
setMethod(
    f = "create_chromatograms",
    signature = c("DBIConnection", "data.frame", "list", "NULL"),
    definition = function(data_obj, metadata, options, features) {
        # TODO: add get_detected_peaks for this
        cd_opts <- options$compound_discoverer

        metadata <- metadata |>
            filter(.data$sample_id %in% options$chromatograms$sample_ids)

        xic_traces_results <- get_xic_traces_from_compounds(
            data_obj,
            cd_opts$compounds_query)

        compounds <- xic_traces_results |>
            group_by(.data$name) |>
            distinct(
                .data$name,
                .data$formula,
                .data$adduct,
                .keep_all = TRUE
            ) |>
            ungroup() |>
            select(name, formula, adduct, mz) |>
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
                    inner_join(compound_data, by = cols) |>
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
                        feature_metadata_id = compound_data$index
                    )

                    additional_metadata_list[[j]] <- data.frame(
                        feature_metadata_id = compound_data$index,
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
                feature_metadata = do.call(rbind, additional_metadata_list),
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

        list(
            chromatograms = do.call(
                rbind,
                lapply(results, `[[`, "chromatograms")),
            mass_traces = data.frame(
                rt = numeric(),
                mz = numeric(),
                metadata_index = numeric(),
                feature_metadata_id = numeric()
            ),
            feature_metadata = do.call(
                rbind,
                lapply(results, `[[`, "feature_metadata")),
            detected_peaks = do.call(
                rbind,
                lapply(results, `[[`, "detected_peaks"))
        )
    }
)

#' @rdname create_chromatograms
setMethod(
    f = "create_chromatograms",
    signature = c("XcmsRawList", "data.frame", "list", "NULL"),
    definition = function(data_obj, metadata, options, features) {
        metadata <- metadata |>
            filter(.data$sample_id %in% options$chromatograms$sample_ids)

        raw_data <- xcmsraw_to_readers(data_obj)

        process_sample <- function(i) {
            sample_metadata <- metadata[i, ]
            raw_obj <- raw_data[[sample_metadata$sample_path]]

            data <- create_bpc_tic(raw_obj, options$chromatograms$aggregation_fun)
            chroms <- data$chromatograms

            data.frame(
                rt = chroms$rt,
                intensity = chroms$intensity,
                metadata_index = sample_metadata$sample_index,
                feature_metadata_id = NA_real_
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

        list(
            chromatograms = do.call(rbind, chromatograms_list),
            mass_traces = data.frame(
                rt = numeric(),
                mz = numeric(),
                metadata_index = numeric(),
                feature_metadata_id = numeric()
            ),
            feature_metadata = data.frame(feature_metadata_id = numeric(), metadata_index = numeric()),
            detected_peaks = data.frame()
        )
    }
)

#' @rdname create_chromatograms
setMethod(
    f = "create_chromatograms",
    signature = c("XcmsRawList", "data.frame", "list", "ANY"),
    definition = function(data_obj, metadata, options, features) {
        metadata <- metadata |>
            filter(.data$sample_id %in% options$chromatograms$sample_ids)

        raw_data <- xcmsraw_to_readers(data_obj)

        process_sample <- function(i) {
            sample_metadata <- metadata[i, ]
            raw_obj <- raw_data[[sample_metadata$sample_path]]

            hdr <- ms_header(raw_obj)
            full_rt_range <- range(hdr$retentionTime)

            feats <- get_features(
                options,
                sample_metadata,
                full_rt_range = full_rt_range
            )

            chromatograms_list <- list()
            mass_traces_list <- list()
            additional_metadata_list <- list()

            n_features <- nrow(options$chromatograms$features)

            for (j in seq_len(length(feats))) {
                feature_data <- feats[[j]]

                data <- create_chromatogram(
                    raw_obj,
                    mz_range = feature_data$mzr,
                    rt_range = feature_data$rtr,
                    fill_gaps = options$chromatograms$fill_gaps
                )

                feature_metadata_id <- (i - 1) * n_features + j

                additional_metadata_list[[j]] <- data.frame(
                    feature_metadata_id = feature_metadata_id,
                    metadata_index = sample_metadata$sample_index,
                    feature_id = feature_data$feature_id
                )

                chromatograms_list[[j]] <- data.frame(
                    rt = data$chromatograms$rt,
                    intensity = data$chromatograms$intensity,
                    metadata_index = sample_metadata$sample_index,
                    feature_metadata_id = feature_metadata_id
                )

                if (nrow(data$mass_traces) > 0) {
                    mass_traces_list[[j]] <- data.frame(
                        rt = data$mass_traces$rt,
                        mz = data$mass_traces$mz,
                        metadata_index = sample_metadata$sample_index,
                        feature_metadata_id = feature_metadata_id
                    )
                }
            }

            if (length(mass_traces_list) == 0) {
                mass_traces <- data.frame(
                    rt = numeric(),
                    mz = numeric(),
                    metadata_index = numeric(),
                    feature_metadata_id = numeric()
                )
            } else {
                mass_traces <- do.call(rbind, mass_traces_list)
            }

            list(
                chromatograms = do.call(rbind, chromatograms_list),
                mass_traces = mass_traces,
                feature_metadata = do.call(rbind, additional_metadata_list)
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

        list(
            chromatograms = do.call(
                rbind,
                lapply(results, `[[`, "chromatograms")),
            mass_traces = do.call(
                rbind,
                lapply(results, `[[`, "mass_traces")),
            feature_metadata = do.call(
                rbind,
                lapply(results, `[[`, "feature_metadata")),
            detected_peaks = data.frame()
        )
    }
)

#' @rdname create_chromatograms
setMethod(
    f = "create_chromatograms",
    signature = c("MChromatograms", "data.frame", "list", "NULL"),
    definition = function(data_obj, metadata, options, features) {
        metadata <- metadata |>
            filter(.data$sample_id %in% options$chromatograms$sample_ids)

        chromatograms <- data.frame()

        for (j in seq_len(ncol(data_obj))) {
            sample_metadata <- metadata |>
                filter(.data$sample_index == j)

            if (nrow(sample_metadata) == 0) next

            for (i in seq_len(nrow(data_obj))) {
                chrom <- data_obj[i, j]

                chromatograms <- rbind(chromatograms, data.frame(
                    rt = MSnbase::rtime(chrom),
                    intensity = MSnbase::intensity(chrom),
                    metadata_index = sample_metadata$sample_index,
                    feature_metadata_id = NA_real_
                ))
            }
        }

        list(
            chromatograms = chromatograms,
            mass_traces = data.frame(
                rt = numeric(),
                mz = numeric(),
                metadata_index = numeric(),
                feature_metadata_id = numeric()
            ),
            feature_metadata = data.frame(feature_metadata_id = numeric(), metadata_index = numeric()),
            detected_peaks = data.frame()
        )
    }
)

#' @rdname create_chromatograms
setMethod(
    f = "create_chromatograms",
    signature = c("XChromatogram", "data.frame", "list", "NULL"),
    definition = function(data_obj, metadata, options, features) {
        chromatograms <- data.frame(
            rt = MSnbase::rtime(data_obj),
            intensity = MSnbase::intensity(data_obj),
            metadata_index = 1L,
            feature_metadata_id = NA_real_
        )

        detected_peaks_raw <- get_detected_peaks(data_obj)
        if (!is.null(detected_peaks_raw)) {
            detected_peaks <- detected_peaks_raw |>
                left_join(metadata, by = "sample_index")
        } else {
            detected_peaks <- data.frame()
        }

        list(
            chromatograms = chromatograms,
            mass_traces = data.frame(
                rt = numeric(),
                mz = numeric(),
                metadata_index = numeric(),
                feature_metadata_id = numeric()
            ),
            feature_metadata = data.frame(feature_metadata_id = numeric(), metadata_index = numeric()),
            detected_peaks = detected_peaks
        )
    }
)

#' @rdname create_chromatograms
setMethod(
    f = "create_chromatograms",
    signature = c("XChromatograms", "data.frame", "list", "NULL"),
    definition = function(data_obj, metadata, options, features) {
        metadata <- metadata |>
            filter(.data$sample_id %in% options$chromatograms$sample_ids)

        n_features <- nrow(data_obj)

        # Pre-compute feature_id per row: mz and RT range are the same across
        # all samples for a given row, so derive them from column 1.
        row_feature_info <- lapply(seq_len(n_features), function(i) {
            mzr <- MSnbase::mz(data_obj[i, 1L])[[1L]]
            rts <- MSnbase::rtime(data_obj[i, 1L])
            mz_center <- mean(mzr)
            rt_center <- if (length(rts) > 0L) mean(range(rts)) else NA_real_
            feature_id <- if (!is.na(rt_center)) {
                sprintf("M%dT%d", round(mz_center), round(rt_center))
            } else {
                sprintf("M%d", round(mz_center))
            }
            list(feature_id = feature_id)
        })

        chromatograms_list <- list()
        feature_metadata_list <- list()
        counter <- 0L

        for (s in seq_len(nrow(metadata))) {
            sample_metadata <- metadata[s, ]
            j <- sample_metadata$sample_index

            for (i in seq_len(n_features)) {
                counter <- counter + 1L
                chrom <- data_obj[i, j]

                chromatograms_list[[counter]] <- data.frame(
                    rt = MSnbase::rtime(chrom),
                    intensity = MSnbase::intensity(chrom),
                    metadata_index = sample_metadata$sample_index,
                    feature_metadata_id = counter
                )

                feature_metadata_list[[counter]] <- data.frame(
                    feature_metadata_id = counter,
                    metadata_index = sample_metadata$sample_index,
                    feature_id = row_feature_info[[i]]$feature_id
                )
            }
        }

        if (length(chromatograms_list) == 0L) {
            chromatograms <- data.frame(
                rt = numeric(),
                intensity = numeric(),
                metadata_index = numeric(),
                feature_metadata_id = numeric()
            )
            feature_metadata <- data.frame(
                feature_metadata_id = numeric(),
                metadata_index = numeric(),
                feature_id = character()
            )
        } else {
            chromatograms <- do.call(rbind, chromatograms_list)
            feature_metadata <- do.call(rbind, feature_metadata_list)
        }

        detected_peaks_raw <- get_detected_peaks(data_obj)
        if (!is.null(detected_peaks_raw)) {
            detected_peaks <- detected_peaks_raw |>
                left_join(metadata, by = "sample_index")
        } else {
            detected_peaks <- data.frame()
        }

        list(
            chromatograms = chromatograms,
            mass_traces = data.frame(
                rt = numeric(),
                mz = numeric(),
                metadata_index = numeric(),
                feature_metadata_id = numeric()
            ),
            feature_metadata = feature_metadata,
            detected_peaks = detected_peaks
        )
    }
)

#' @rdname create_chromatograms
setMethod(
    f = "create_chromatograms",
    signature = c("ANY", "data.frame", "list", "NULL"),
    definition = function(data_obj, metadata, options, features) {
        metadata <- metadata |>
            filter(.data$sample_id %in% options$chromatograms$sample_ids)
        raw_data <- io_get_raw_data(metadata$sample_path)

        process_sample <- function(i) {
            sample_metadata <- metadata[i, ]
            raw_obj <- raw_data[[sample_metadata$sample_path]]

            if (options$chromatograms$rt_type %in% c("corrected", "both")) {
                if (!is_xcms_data(data_obj)) {
                    stop("The data object should be XCMSnExp or MsExperiment to plot the RT adjusted chromatograms")
                }

                scan_indices <- which(
                    xcms::fromFile(data_obj) == sample_metadata$sample_index
                )
                rt_adjusted <- xcms::adjustedRtime(data_obj)[scan_indices]
            } else {
                rt_adjusted <- NULL
            }

            data <- create_bpc_tic(
                raw_obj,
                options$chromatograms$aggregation_fun
            )
            chroms <- data$chromatograms

            if (options$chromatograms$rt_type == "both") {
                chroms_raw <- chroms |> mutate(rt_type = "RT uncorrected")
                chroms_adj <- chroms |> mutate(rt = rt_adjusted, rt_type = "RT corrected")
                chroms <- rbind(chroms_raw, chroms_adj)
            } else if (options$chromatograms$rt_type == "corrected") {
                chroms <- chroms |> mutate(rt = rt_adjusted, rt_type = "RT corrected")
            } else {
                chroms <- chroms |> mutate(rt_type = NA_character_)
            }

            data.frame(
                rt = chroms$rt,
                intensity = chroms$intensity,
                rt_type = chroms$rt_type,
                metadata_index = sample_metadata$sample_index,
                feature_metadata_id = NA_real_
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

        list(
            chromatograms = chromatograms,
            mass_traces = data.frame(
                rt = numeric(),
                mz = numeric(),
                metadata_index = numeric(),
                feature_metadata_id = numeric()
            ),
            feature_metadata = data.frame(feature_metadata_id = numeric(), metadata_index = numeric()),
            detected_peaks = data.frame()
        )
    }
)

#' @rdname create_chromatograms
setMethod(
    f = "create_chromatograms",
    signature = c("ANY", "data.frame", "list", "character"),
    definition = function(data_obj, metadata, options, features) {
        if (!is_xcms_data(data_obj)) {
            stop("To use feature IDs from the grouped peaks you need to provide an xcms object.")
        }

        metadata <- metadata |>
            filter(.data$sample_id %in% options$chromatograms$sample_ids)
        raw_data <- io_get_raw_data(metadata$sample_path)
        adjusted_rts <- get_adjusted_rts(data_obj)
        all_detected_peaks <- get_detected_peaks(data_obj)
        grouped_peaks <- get_grouped_peaks(data_obj) |>
            filter(.data$name %in% options$chromatograms$features)
        detected_peaks <- data.frame()

        chromatograms <- data.frame()
        mass_traces <- data.frame(
            rt = numeric(),
            mz = numeric(),
            metadata_index = numeric(),
            feature_metadata_id = numeric()
        )
        feature_metadata <- data.frame(feature_metadata_id = numeric(), metadata_index = numeric())
        feature_metadata_counter <- 0L

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

                feature_metadata_counter <- feature_metadata_counter + 1L
                feature_metadata <- rbind(feature_metadata, data.frame(
                    feature_metadata_id = feature_metadata_counter,
                    metadata_index = sample_metadata$sample_index,
                    feature_id = feature$name
                ))

                chromatograms <- rbind(chromatograms, data.frame(
                    rt = data$chromatograms$rt,
                    intensity = data$chromatograms$intensity,
                    metadata_index = sample_metadata$sample_index,
                    feature_metadata_id = feature_metadata_counter
                ))

                if (nrow(data$mass_traces) > 0) {
                    mass_traces <- rbind(mass_traces, data.frame(
                        rt = data$mass_traces$rt,
                        mz = data$mass_traces$mz,
                        metadata_index = sample_metadata$sample_index,
                        feature_metadata_id = feature_metadata_counter
                    ))
                }
            }
        }

        io_close_raw_data(raw_data)

        list(
            chromatograms = chromatograms,
            mass_traces = mass_traces,
            feature_metadata = feature_metadata,
            detected_peaks = detected_peaks
        )
    }
)

#' @rdname create_chromatograms
setMethod(
    f = "create_chromatograms",
    signature = c("ANY", "data.frame", "list", "ANY"),
    definition = function(data_obj, metadata, options, features) {
        metadata <- metadata |>
            filter(.data$sample_id %in% options$chromatograms$sample_ids)

        raw_data <- io_get_raw_data(metadata$sample_path)
        all_detected_peaks <- get_detected_peaks(data_obj)
        adjusted_rts <- get_adjusted_rts(data_obj)

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

            feats <- get_features(
                options,
                sample_metadata,
                full_rt_range = full_rt_range
            )

            for (j in seq_len(length(feats))) {
                feature_data <- feats[[j]]

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
                feature_metadata_id <- (i - 1) * n_features + j

                additional_metadata_list[[j]] <- data.frame(
                    feature_metadata_id = feature_metadata_id,
                    metadata_index = sample_metadata$sample_index,
                    feature_id = feature_data$feature_id
                )

                chromatograms_list[[j]] <- data.frame(
                    rt = data$chromatograms$rt,
                    intensity = data$chromatograms$intensity,
                    metadata_index = sample_metadata$sample_index,
                    feature_metadata_id = feature_metadata_id
                )

                if (nrow(data$mass_traces) > 0) {
                    mass_traces_list[[j]] <- data.frame(
                        rt = data$mass_traces$rt,
                        mz = data$mass_traces$mz,
                        metadata_index = sample_metadata$sample_index,
                        feature_metadata_id = feature_metadata_id
                    )
                }
            }

            if (length(mass_traces_list) == 0) {
                mass_traces <- data.frame(
                    rt = numeric(),
                    mz = numeric(),
                    metadata_index = numeric(),
                    feature_metadata_id = numeric()
                )
            } else {
                mass_traces <- do.call(rbind, mass_traces_list)
            }

            list(
                chromatograms = do.call(rbind, chromatograms_list),
                mass_traces = mass_traces,
                feature_metadata = do.call(rbind, additional_metadata_list),
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

        list(
            chromatograms = do.call(
                rbind,
                lapply(results, `[[`, "chromatograms")),
            mass_traces = do.call(
                rbind,
                lapply(results, `[[`, "mass_traces")),
            feature_metadata = do.call(
                rbind,
                lapply(results, `[[`, "feature_metadata")),
            detected_peaks = do.call(
                rbind,
                lapply(results, `[[`, "detected_peaks"))
        )
    }
)
