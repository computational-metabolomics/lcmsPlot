#' Create the spectra of a data container
#'
#' Every spectra mode is reduced to one table of target scans (by scan index,
#' or by the retention time closest to a target), from which the spectra are
#' taken in a single subset of the raw `Spectra` and held in memory. Spectra
#' read at detected peaks carry the peak's `feature_id`, so they line up with
#' the chromatograms the peaks were matched to.
#'
#' Intensities are scaled to a maximum of 100 per spectrum. When a spectral
#' library is given, the best matching library spectrum of each spectrum is
#' added with `reference = TRUE`.
#'
#' @param obj An `lcmsPlotDataContainer`.
#' @param options A `list` representing the plot object's options.
#' @return The `lcmsPlotDataContainer` with its `spectra` slot set: a
#' `Spectra` object with spectra variables `sample_index`, `feature_id`, and
#' `reference`.
#' @keywords internal
create_spectra <- function(obj, options) {
    opt <- options$spectra
    is_standalone <- !options$chromatograms$show & opt$show
    samples <- .select_samples(obj, opt$sample_ids)
    obj <- ensure_spectra(obj, options$parallel_param)

    targets <- spectra_targets(obj, samples, opt, is_standalone)
    sps <- select_target_spectra(obj, targets, opt$ms_level)

    if (!is.null(opt$spectral_match_db) && length(sps) > 0) {
        sps <- c(sps, match_library_spectra(
            sps,
            read_spectral_library(opt$spectral_match_db),
            opt$match_target_index %||% 1L))
    }

    obj@spectra <- scale_spectra(sps)
    obj
}

#' Determine the scans to read the spectra from
#'
#' @param obj An `lcmsPlotDataContainer`.
#' @param samples A `data.frame` of the samples to read spectra from.
#' @param opt The `spectra` options.
#' @param is_standalone A `logical` indicating whether the spectra are plotted
#' without chromatograms.
#' @return A `tibble` with columns `sample_index`, `feature_id`, and either
#' `scan_index` (the position of the scan in its sample) or `target_rt`.
#' @keywords internal
spectra_targets <- function(obj, samples, opt, is_standalone) {
    if (!is.null(opt$scan_index)) {
        scan_index <- if (is.character(opt$scan_index)) {
            samples[[opt$scan_index]]
        } else {
            rep_len(opt$scan_index, nrow(samples))
        }
        return(tibble(
            sample_index = samples$sample_index,
            feature_id = NA_character_,
            scan_index = as.integer(scan_index)))
    }

    if (opt$mode == "closest") {
        if (is.null(opt$rt)) {
            stop("lp_spectra: mode = \"closest\" requires 'rt'.")
        }
        # Different RTs may be supplied per sample ID.
        target_rt <- if (!is.null(names(opt$rt))) {
            unname(opt$rt[samples$sample_id])
        } else {
            rep_len(opt$rt, nrow(samples))
        }
        return(tibble(
            sample_index = samples$sample_index,
            feature_id = NA_character_,
            target_rt = target_rt))
    }

    # With chromatograms, the peaks they selected (and matched to features);
    # standalone, every detected peak of the requested samples.
    peaks <- obj@detected_peaks
    if (is_standalone && nrow(peaks) == 0) {
        peaks <- .source_peaks(obj, samples)
    }
    if (nrow(peaks) == 0) {
        return(tibble(
            sample_index = integer(),
            feature_id = character(),
            target_rt = numeric()))
    }

    peaks <- peaks |>
        filter(.data$sample_index %in% samples$sample_index)
    feature_id <- if ("feature_id" %in% colnames(peaks)) {
        as.character(peaks$feature_id)
    } else {
        rep(NA_character_, nrow(peaks))
    }

    if (opt$mode == "closest_apex") {
        return(tibble(
            sample_index = peaks$sample_index,
            feature_id = feature_id,
            target_rt = peaks$rt))
    }

    # across_peak: one spectrum every `interval` seconds across each peak.
    rts <- lapply(seq_len(nrow(peaks)), function(i) {
        seq(peaks$rtmin[i], peaks$rtmax[i], by = opt$interval)
    })
    n_rts <- lengths(rts)
    tibble(
        sample_index = rep(peaks$sample_index, n_rts),
        feature_id = rep(feature_id, n_rts),
        target_rt = unlist(rts))
}

#' Read the spectra of the target scans
#'
#' @param obj An `lcmsPlotDataContainer` whose spectra are loaded.
#' @param targets A `data.frame` as returned by [spectra_targets()].
#' @param ms_level The MS level of the scans matched by retention time.
#' @return An in-memory `Spectra` object, one spectrum per target, with
#' spectra variables `sample_index`, `feature_id`, and `reference`.
#' @keywords internal
select_target_spectra <- function(obj, targets, ms_level) {
    if (nrow(targets) == 0) {
        return(Spectra::Spectra())
    }

    by_sample <- split(
        seq_len(nrow(targets)),
        factor(targets$sample_index, levels = unique(targets$sample_index)))

    parts <- lapply(names(by_sample), function(sample_key) {
        rows <- by_sample[[sample_key]]
        target <- targets[rows, , drop = FALSE]

        if ("scan_index" %in% colnames(target)) {
            sps <- raw_spectra(obj, as.integer(sample_key))
            idx <- target$scan_index
        } else {
            sps <- raw_spectra(obj, as.integer(sample_key), ms_level)
            if (length(sps) == 0) {
                stop("No scans found with the specified MS level.")
            }
            idx <- .closest_scans(Spectra::rtime(sps), target$target_rt)
        }

        selected <- sps[idx]
        selected$feature_id <- target$feature_id
        selected
    })

    sps <- Spectra::setBackend(
        do.call(c, unname(parts)), Spectra::MsBackendMemory())
    sps$reference <- FALSE
    sps
}

# The index of the scan closest to each target retention time.
.closest_scans <- function(scan_rt, target_rt) {
    scan_order <- order(scan_rt)
    target_order <- order(target_rt)
    idx <- integer(length(target_rt))
    idx[target_order] <- scan_order[MsCoreUtils::closest(
        target_rt[target_order],
        scan_rt[scan_order],
        tolerance = Inf,
        duplicates = "keep")]
    idx
}

read_spectral_library <- function(paths) {
    source <- if (length(paths) == 1 && all(endsWith(paths, ".msp"))) {
        MsBackendMsp::MsBackendMsp()
    } else {
        Spectra::MsBackendMzR()
    }
    Spectra::Spectra(paths, source = source)
}

#' Find the library spectra matching each query spectrum
#'
#' @param query A `Spectra` object with spectra variables `sample_index` and
#' `feature_id`.
#' @param library A `Spectra` object holding the spectral library.
#' @param match_target_index An `integer` selecting which match to keep: `1`
#' for the best match, `2` for the second best, and so on.
#' @return An in-memory `Spectra` with the matching library spectra, carrying
#' the `rtime`, `sample_index`, and `feature_id` of their query spectrum and
#' `reference = TRUE`.
#' @keywords internal
match_library_spectra <- function(query, library, match_target_index = 1L) {
    similarities <- Spectra::compareSpectra(query, library)
    if (is.null(dim(similarities))) {
        similarities <- matrix(similarities, nrow = length(query))
    }

    hits <- lapply(seq_len(nrow(similarities)), function(i) {
        values <- sort(unique(similarities[i, ]), decreasing = TRUE)
        if (length(values) < match_target_index) {
            return(integer())
        }
        which(similarities[i, ] == values[match_target_index])
    })
    query_index <- rep(seq_along(hits), lengths(hits))

    matches <- Spectra::setBackend(
        library[unlist(hits)], Spectra::MsBackendMemory())
    Spectra::rtime(matches) <- Spectra::rtime(query)[query_index]
    matches$sample_index <- query$sample_index[query_index]
    matches$feature_id <- query$feature_id[query_index]
    matches$reference <- TRUE
    matches
}

# Scale every spectrum's intensities to a maximum of 100.
scale_spectra <- function(sps) {
    if (length(sps) == 0) {
        return(sps)
    }
    sps <- Spectra::addProcessing(sps, function(x, ...) {
        max_intensity <- suppressWarnings(max(x[, 2L], na.rm = TRUE))
        if (is.finite(max_intensity) && max_intensity > 0) {
            x[, 2L] <- 100 * x[, 2L] / max_intensity
        }
        x
    })
    Spectra::applyProcessing(sps)
}
