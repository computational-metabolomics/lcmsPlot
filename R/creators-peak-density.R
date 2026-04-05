#' Create an instance of class `lcmsPlotDataContainer` from a peak density
#' dataset.
#'
#' This function creates an `lcmsPlotDataContainer` object from
#' a dataset that contains peak density data for each feature (i.e.,
#' the kernel density of peak apex RTs across samples for each m/z bin),
#' mirroring the algorithm used by `xcms::plotChromPeakDensity()`.
#'
#' @param obj An instance of class `lcmsPlotDataContainer`.
#' @param options A `list` representing the plot object's options.
#' @return An instance of class `lcmsPlotDataContainer` with the created
#' peak density dataset. The `peak_density` slot is a `tibble` with columns
#' `rt`, `density`, `rtmin`, `rtmax`, `data_type` (`"density"` or `"rect"`),
#' `mzmin`, `mzmax`, `metadata_index`, and `feature_metadata_id`.
#' @keywords internal
setGeneric(
    "create_peak_density",
    function(obj, options) standardGeneric("create_peak_density")
)

# Descend from the local maximum to find feature boundaries.
# Mirrors xcms::descendMin used in the peak density grouping algorithm.
.descend_min <- function(y, idx) {
    n <- length(y)
    left <- idx
    while (left > 1L && y[left - 1L] <= y[left]) left <- left - 1L
    right <- idx
    while (right < n && y[right + 1L] <= y[right]) right <- right + 1L
    c(left, right)
}

#' @rdname create_peak_density
setMethod(
    f = "create_peak_density",
    signature = c("lcmsPlotDataContainer", "list"),
    definition = function(obj, options) {
        all_peaks <- as_tibble(as.data.frame(xcms::chromPeaks(obj@data_obj))) |>
            dplyr::rename(sample_index = sample)

        # Populate detected_peaks from data_obj if not already done
        if (nrow(obj@detected_peaks) == 0) {
            obj@detected_peaks <- all_peaks |>
                left_join(obj@metadata, by = "sample_index")
        }

        opts <- options$peak_density
        features <- opts$features
        bw <- opts$bw
        min_fraction <- opts$min_fraction
        min_samples <- opts$min_samples
        max_features <- opts$max_features
        n_samples <- nrow(obj@metadata)

        sample_groups <- if (is.null(opts$sample_groups)) {
            rep(1L, n_samples)
        } else {
            opts$sample_groups
        }
        sample_groups_table <- table(sample_groups)

        # Density range computed over ALL peaks, matching xcms exactly
        all_rt <- all_peaks$rt
        full_rt_range <- range(all_rt)
        dens_from <- full_rt_range[1] - 3 * bw
        dens_to   <- full_rt_range[2] + 3 * bw
        densN <- max(512L,
                     2L * 2L^ceiling(log2(diff(full_rt_range) / (bw / 2))))

        simulate <- !is.null(min_fraction)

        results <- lapply(seq_len(nrow(features)), function(i) {
            feat <- features[i, ]

            peaks_in <- all_peaks |>
                filter(
                    .data$mz >= feat[["mzmin"]],
                    .data$mz <= feat[["mzmax"]]
                )

            has_rt <- !is.na(feat[["rtmin"]]) && !is.na(feat[["rtmax"]])
            if (has_rt) {
                peaks_in <- peaks_in |>
                    filter(
                        .data$rt >= feat[["rtmin"]],
                        .data$rt <= feat[["rtmax"]]
                    )
            }

            if (nrow(peaks_in) == 0) {
                return(NULL)
            }

            dens <- density(
                peaks_in$rt,
                bw = bw,
                from = dens_from,
                to = dens_to,
                n = densN
            )

            # Clip the density curve to the feature's RT window so the plot
            # does not display the tails beyond the specified range.
            clip_from <- if (has_rt) feat[["rtmin"]] else min(peaks_in$rt)
            clip_to <- if (has_rt) feat[["rtmax"]] else max(peaks_in$rt)
            keep <- dens$x >= clip_from & dens$x <= clip_to

            density_rows <- tibble(
                rt = dens$x[keep],
                density = dens$y[keep],
                rtmin = NA_real_,
                rtmax = NA_real_,
                data_type = "density",
                mzmin = feat[["mzmin"]],
                mzmax = feat[["mzmax"]],
                metadata_index = NA_real_,
                feature_metadata_id = as.integer(i)
            )

            rect_rows <- if (simulate) {
                dens_y <- dens$y
                dens_max <- max(dens_y)
                snum <- 0L
                rects <- list()

                while (dens_y[max_y <- which.max(dens_y)] > dens_max / 20 &&
                       snum < max_features) {
                    feat_range <- .descend_min(dens_y, max_y)
                    dens_y[feat_range[1]:feat_range[2]] <- 0

                    feat_idx <- which(
                        peaks_in$rt >= dens$x[feat_range[1]] &
                        peaks_in$rt <= dens$x[feat_range[2]]
                    )

                    if (length(feat_idx) == 0) next

                    tt <- table(
                        sample_groups[unique(peaks_in$sample_index[feat_idx])]
                    )
                    passes <- any(
                        tt / sample_groups_table[names(tt)] >= min_fraction &
                        tt >= min_samples
                    )
                    if (!passes) next

                    rects[[length(rects) + 1L]] <- tibble(
                        rt = NA_real_,
                        density = NA_real_,
                        rtmin = min(peaks_in$rt[feat_idx]),
                        rtmax = max(peaks_in$rt[feat_idx]),
                        data_type = "rect",
                        mzmin = feat[["mzmin"]],
                        mzmax = feat[["mzmax"]],
                        metadata_index = NA_real_,
                        feature_metadata_id = as.integer(i)
                    )
                    snum <- snum + 1L
                }

                bind_rows(rects)
            } else {
                tibble()
            }

            bind_rows(density_rows, rect_rows)
        })

        obj@peak_density <- bind_rows(results)

        validObject(obj)

        return(obj)
    }
)
