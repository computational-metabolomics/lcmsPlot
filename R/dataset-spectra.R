#' Create the spectra of a data container
#'
#' Every spectra mode is reduced to one table of target scans (by scan index,
#' or by the retention time closest to a target), from which the spectra are
#' taken in a single subset of the raw `Spectra` and held in memory. Spectra
#' read at detected peaks carry the peak's `feature_id`, so they line up with
#' the chromatograms the peaks were matched to. When the data has adjusted
#' retention times, as the detected peaks do, these spectra are selected and
#' reported on them; spectra read by scan index or at a given RT are not.
#'
#' Intensities are scaled to a maximum of 100 per spectrum. When a spectral
#' library is given, the best matching library spectrum of each spectrum is
#' added with `reference = TRUE`. MS2 spectra carry their precursor m/z in a
#' `precursor_mz` spectra variable.
#'
#' @param obj An `lcmsPlotDataContainer`.
#' @param options A `list` representing the plot object's options.
#' @return The `lcmsPlotDataContainer` with its `spectra` slot set: a
#' `Spectra` object with spectra variables `sample_index`, `feature_id`,
#' `reference`, and, for MS2 spectra, `precursor_mz`.
#' @keywords internal
create_spectra <- function(obj, options) {
    opt <- options$spectra
    is_standalone <- !options$chromatograms$show & opt$show
    samples <- .select_samples(obj, opt$sample_ids)
    obj <- ensure_spectra(obj, options$parallel_param)
    by_scan_index <- !is.null(opt$scan_index)
    by_precursor <- !by_scan_index && opt$mode == "precursor"
    # Detected peaks are on the adjusted retention times, so the scans read
    # at peaks are matched and reported on them too.
    at_peaks <- !by_scan_index && opt$mode != "closest"

    targets <- spectra_targets(obj, samples, opt, is_standalone)
    sps <- select_target_spectra(
        obj, targets, opt$ms_level,
        adjusted = at_peaks && has_adjusted_rtime(obj))

    if (by_precursor && opt$select == "consensus" && length(sps) > 0) {
        sps <- consensus_spectra(
            sps, targets$group, opt$consensus_ppm, opt$consensus_min_prop)
    }

    # Some readers report a precursor m/z of 0 for MS1 scans.
    is_msn <- Spectra::msLevel(sps) > 1L
    if (any(is_msn)) {
        sps$precursor_mz <- ifelse(
            is_msn, Spectra::precursorMz(sps), NA_real_)
    }

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
#' Precursor targets add a `group` column, one value per feature and sample.
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

    if (opt$mode == "precursor") {
        return(precursor_targets(obj, samples, opt, is_standalone))
    }

    peaks <- .target_peaks(obj, samples, is_standalone)
    if (nrow(peaks) == 0) {
        return(tibble(
            sample_index = integer(),
            feature_id = character(),
            target_rt = numeric()))
    }

    feature_id <- .peak_feature_ids(peaks)

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

# With chromatograms, the peaks they selected (and matched to features);
# standalone, every detected peak of the requested samples.
.target_peaks <- function(obj, samples, is_standalone) {
    peaks <- obj@detected_peaks
    if (is_standalone && nrow(peaks) == 0) {
        peaks <- .source_peaks(obj, samples)
    }
    if (nrow(peaks) == 0) {
        return(peaks)
    }
    peaks |> filter(.data$sample_index %in% samples$sample_index)
}

.peak_feature_ids <- function(peaks) {
    if ("feature_id" %in% colnames(peaks)) {
        as.character(peaks$feature_id)
    } else {
        rep(NA_character_, nrow(peaks))
    }
}

#' Find the MS2 scans fragmenting the precursor of each feature
#'
#' A scan matches a feature in a sample when its precursor m/z lies within
#' `opt$ppm` of the m/z of the feature's detected peak and its retention time
#' within the peak's `[rtmin, rtmax]`. A feature with no detected peak in a
#' sample uses the extraction window of its chromatogram instead. Only the
#' scan headers are read. Features with no matching scan are skipped with a
#' single warning listing them.
#'
#' @param obj An `lcmsPlotDataContainer` whose spectra are loaded.
#' @param samples A `data.frame` of the samples to read spectra from.
#' @param opt The `spectra` options; uses `ppm` and `select`.
#' @param is_standalone A `logical` indicating whether the spectra are plotted
#' without chromatograms.
#' @return A `tibble` with columns `sample_index`, `feature_id`, `scan_index`
#' (the position of the scan in its sample), and `group`, ordered by sample.
#' @keywords internal
precursor_targets <- function(obj, samples, opt, is_standalone) {
    windows <- .precursor_windows(obj, samples, opt$ppm, is_standalone)
    adjusted <- has_adjusted_rtime(obj)

    by_sample <- split(windows, factor(
        windows$sample_index, levels = samples$sample_index))
    targets <- bind_rows(lapply(by_sample, function(w) {
        if (nrow(w) == 0) {
            return(NULL)
        }
        sps <- select_rtime(raw_spectra(obj, w$sample_index[1]), adjusted)
        .select_precursor_scans(sps, w, opt$select)
    }))
    if (nrow(targets) == 0) {
        targets <- tibble(
            sample_index = integer(),
            feature_id = character(),
            scan_index = integer(),
            group = character())
    }

    missing <- windows |>
        distinct(.data$group, .keep_all = TRUE) |>
        filter(!.data$group %in% targets$group) |>
        arrange(match(.data$sample_index, samples$sample_index))
    if (nrow(missing) > 0) {
        sample_ids <- samples$sample_id[
            match(missing$sample_index, samples$sample_index)]
        warning(
            "lp_spectra: no MS2 scan matches the precursor of ",
            nrow(missing), " feature(s), skipped: ",
            paste0(missing$feature_id, " (", sample_ids, ")", collapse = ", "),
            call. = FALSE)
    }
    targets
}

# The precursor windows of every feature and sample: one per detected peak,
# plus the chromatogram's extraction window for the features with no peak in
# a sample.
.precursor_windows <- function(obj, samples, ppm, is_standalone) {
    peaks <- .target_peaks(obj, samples, is_standalone)
    peak_windows <- if (nrow(peaks) > 0) {
        tibble(
            sample_index = as.integer(peaks$sample_index),
            feature_id = .peak_feature_ids(peaks),
            mzmin = peaks$mz - peaks$mz * ppm / 1e6,
            mzmax = peaks$mz + peaks$mz * ppm / 1e6,
            rtmin = peaks$rtmin,
            rtmax = peaks$rtmax,
            apex_rt = peaks$rt)
    } else {
        tibble(
            sample_index = integer(),
            feature_id = character(),
            mzmin = numeric(),
            mzmax = numeric(),
            rtmin = numeric(),
            rtmax = numeric(),
            apex_rt = numeric())
    }
    peak_windows <- peak_windows |> filter(.complete_windows(peak_windows))

    chrom_data <- as.data.frame(Chromatograms::chromData(obj@chromatograms))
    window_cols <- c(
        "sample_index", "feature_id", "mzMin", "mzMax", "rtMin", "rtMax")
    feature_windows <- if (all(window_cols %in% colnames(chrom_data))) {
        tibble(
            sample_index = as.integer(chrom_data$sample_index),
            feature_id = as.character(chrom_data$feature_id),
            mzmin = chrom_data$mzMin,
            mzmax = chrom_data$mzMax,
            rtmin = chrom_data$rtMin,
            rtmax = chrom_data$rtMax,
            apex_rt = (chrom_data$rtMin + chrom_data$rtMax) / 2) |>
            filter(
                .data$sample_index %in% samples$sample_index,
                !is.na(.data$feature_id)) |>
            anti_join(peak_windows, by = c("sample_index", "feature_id"))
    }

    windows <- bind_rows(peak_windows, feature_windows)
    windows <- windows |> filter(.complete_windows(windows))

    # A feature's peaks in one sample are pooled; a peak matched to no feature
    # stands on its own.
    windows$group <- as.character(ifelse(
        is.na(windows$feature_id),
        paste0("peak", seq_len(nrow(windows))),
        paste(windows$sample_index, windows$feature_id)))
    windows
}

.complete_windows <- function(windows) {
    stats::complete.cases(
        windows[, c("mzmin", "mzmax", "rtmin", "rtmax"), drop = FALSE])
}

# The MS2 scans of one sample's spectra matching each group of precursor
# windows, reduced to the scans kept by `select`.
.select_precursor_scans <- function(sps, windows, select) {
    precursor_mz <- Spectra::precursorMz(sps)
    rt <- Spectra::rtime(sps)
    is_candidate <- Spectra::msLevel(sps) == 2L & !is.na(precursor_mz)

    matches <- lapply(seq_len(nrow(windows)), function(k) {
        which(is_candidate &
            precursor_mz >= windows$mzmin[k] &
            precursor_mz <= windows$mzmax[k] &
            rt >= windows$rtmin[k] &
            rt <= windows$rtmax[k])
    })
    by_group <- split(seq_len(nrow(windows)), factor(
        windows$group, levels = unique(windows$group)))

    bind_rows(lapply(by_group, function(rows) {
        idx <- sort(unique(unlist(matches[rows])))
        if (length(idx) == 0) {
            return(NULL)
        }
        idx <- switch(
            select,
            most_intense = idx[which.max(.precursor_intensity(sps[idx]))],
            closest_apex = idx[which.min(vapply(
                rt[idx],
                function(x) min(abs(x - windows$apex_rt[rows])),
                numeric(1)))],
            idx[order(rt[idx])])
        tibble(
            sample_index = windows$sample_index[rows[1]],
            feature_id = windows$feature_id[rows[1]],
            scan_index = idx,
            group = windows$group[rows[1]])
    }))
}

# The intensity of each scan's precursor, or the scans' total ion current when
# the file does not report it for every scan.
.precursor_intensity <- function(sps) {
    intensity <- Spectra::precursorIntensity(sps)
    if (anyNA(intensity) || any(intensity <= 0)) {
        intensity <- .scan_intensity(sps, "sum")
    }
    intensity
}

#' Combine the spectra of each group into one consensus spectrum
#'
#' Keeps the fragments found, within `ppm`, in at least `min_prop` of a
#' group's spectra. The consensus spectrum keeps the spectra variables of the
#' group's first spectrum and the median retention time of its spectra.
#'
#' @param sps A `Spectra` object.
#' @param group A `character` vector with the group of each spectrum.
#' @param ppm A `numeric` m/z tolerance, in ppm, for matching fragments.
#' @param min_prop A `numeric` value in `[0, 1]`: the minimum proportion of a
#' group's spectra a fragment must be found in.
#' @return An in-memory `Spectra` with one spectrum per group, in the order the
#' groups first appear.
#' @keywords internal
consensus_spectra <- function(sps, group, ppm, min_prop) {
    f <- factor(group, levels = unique(group))
    combined <- Spectra::combineSpectra(
        sps,
        f = f,
        FUN = Spectra::combinePeaksData,
        ppm = ppm,
        tolerance = 0,
        peaks = "intersect",
        minProp = min_prop)
    Spectra::rtime(combined) <- as.numeric(
        tapply(Spectra::rtime(sps), f, stats::median))
    combined
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
#' the `rtime`, `sample_index`, `feature_id`, and `precursor_mz` (when set) of
#' their query spectrum and `reference = TRUE`.
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
    if ("precursor_mz" %in% Spectra::spectraVariables(query)) {
        matches$precursor_mz <- query$precursor_mz[query_index]
    }
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
