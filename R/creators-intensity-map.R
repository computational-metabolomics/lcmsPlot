#' Create an instance of class `lcmsPlotDataContainer` from an intensity map
#'
#' This function creates an `lcmsPlotDataContainer` object from
#' an intensity map which is a point-cloud of m/z and RT points.
#'
#' @param obj An instance of class `lcmsPlotDataContainer`.
#' @param options A `list` representing the plot object's options.
#' @return An instance of class `lcmsPlotDataContainer` with the created
#' intensity map, a `tibble` with columns `mz` and `rt`.
#' @keywords internal
setGeneric(
    "create_intensity_map",
    function(obj, options) standardGeneric("create_intensity_map")
)

#' @rdname create_intensity_map
setMethod(
    f = "create_intensity_map",
    signature = c("lcmsPlotDataContainer", "list"),
    definition = function(obj, options) {
        metadata <- obj@metadata |>
            filter(.data$sample_id %in% options$intensity_maps$sample_ids)

        process_sample <- function(i) {
            sample_metadata <- metadata[i, ]
            ms <- open_raw_reader(sample_metadata$sample_path)
            hdr <- ms_scan_index(ms)

            rt_range <- options$intensity_maps$rt_range
            mz_range <- options$intensity_maps$mz_range

            idx <- which(
                hdr$retentionTime >= rt_range[1] &
                    hdr$retentionTime <= rt_range[2])

            # One read for every scan in the window rather than one per scan.
            in_window <- if (length(idx) > 0) {
                .points_in_mz_window(ms_peaks(ms, idx), mz_range)
            } else {
                list()
            }

            ms_close(ms)

            n_points <- vapply(in_window, nrow, integer(1))
            df <- if (sum(n_points) > 0) {
                tibble(
                    rt = rep(hdr$retentionTime[idx], n_points),
                    mz = unlist(
                        lapply(in_window, function(x) x[, 1]),
                        use.names = FALSE),
                    intensity = unlist(
                        lapply(in_window, function(x) x[, 2]),
                        use.names = FALSE)
                )
            } else {
                NULL
            }

            if (!is.null(df) && nrow(df) > 0) {
                if (identical(options$intensity_maps$geom, "point")) {
                    # A scatter of the raw centroids, as xcms::plotRaw() draws:
                    # binning here would imply continuous coverage that the
                    # instrument never measured.
                    binned <- df
                } else {
                    binned <- df |>
                        mutate(
                            rt = bin_coordinate(
                                .data$rt, options$intensity_maps$bin_rt),
                            mz = bin_coordinate(
                                .data$mz, options$intensity_maps$bin_mz)
                        ) |>
                        group_by(.data$rt, .data$mz) |>
                        summarize(
                            intensity = sum(.data$intensity),
                            .groups = "drop"
                        )
                }

                binned |>
                    mutate(
                        metadata_index = sample_metadata$sample_index,
                        feature_metadata_id = sample_metadata$sample_index
                    )
            } else {
                NULL
            }
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

        obj@intensity_maps <- bind_rows(results)

        validObject(obj)

        return(obj)
    }
)
