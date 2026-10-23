#' Create the chromatograms of a data container
#'
#' Dispatches on the type of the container's source and of `features`. Every
#' method reduces its input to a table of extraction windows (one row per
#' chromatogram) and hands it to [extract_chromatograms()], or converts
#' chromatograms the source already holds.
#'
#' @param source The source object of the container (`container@@source`).
#' @param container An `lcmsPlotDataContainer`.
#' @param options A `list` of plot options.
#' @param features `NULL`, a `character` vector of feature IDs (the row names
#' of `xcms::featureDefinitions()`), or a `matrix`/`data.frame` of feature
#' ranges.
#' @return The `lcmsPlotDataContainer` with its `chromatograms` and
#' `detected_peaks` slots set.
#' @keywords internal
setGeneric(
    "create_chromatograms",
    function(source, container, options, features)
        standardGeneric("create_chromatograms")
)

#' @rdname create_chromatograms
setMethod(
    f = "create_chromatograms",
    signature = c(source = "DBIConnection", features = "NULL"),
    definition = function(source, container, options, features) {
        cd_opts <- options$compound_discoverer
        samples <- .select_samples(
            container, options$chromatograms$sample_ids)

        entries <- get_xic_traces_from_compounds(
            source, cd_opts$compounds_query) |>
            as_tibble() |>
            inner_join(
                samples |> select("sample_id", "sample_index"),
                by = "sample_id")

        # One chromatogram per compound and sample, decoded from its stored
        # trace, with the compound annotations on the chromatogram itself.
        traces <- entries |>
            distinct(
                .data$sample_index, .data$name, .data$formula, .data$adduct,
                .keep_all = TRUE) |>
            arrange(match(.data$sample_index, samples$sample_index))

        peaks <- lapply(seq_len(nrow(traces)), function(i) {
            trace <- parse_trace(traces$trace[[i]])
            if (is.null(trace)) {
                return(.empty_peaks())
            }
            keep <- trace$rt >= traces$rtmin[i] - cd_opts$rt_extend &
                trace$rt <= traces$rtmax[i] + cd_opts$rt_extend
            data.frame(rtime = trace$rt[keep], intensity = trace$intensity[keep])
        })

        chrom_data <- traces |>
            transmute(
                sample_index = .data$sample_index,
                feature_id = .data$name,
                name = .data$name,
                formula = .data$formula,
                adduct = .data$adduct,
                mz = .data$mz,
                rtMin = .data$rtmin - cd_opts$rt_extend,
                rtMax = .data$rtmax + cd_opts$rt_extend
            )

        container@chromatograms <- as_memory_chromatograms(chrom_data, peaks)
        container@detected_peaks <- normalise_peaks(
            entries |>
                mutate(feature_id = .data$name) |>
                select(
                    "name", "feature_id", "mz", "rt", "rtmin", "rtmax",
                    "into", "maxo", "sample_index"))
        container
    }
)

#' @rdname create_chromatograms
setMethod(
    f = "create_chromatograms",
    signature = c(source = "CompoundDiscovererNodeSource", features = "NULL"),
    definition = function(source, container, options, features) {
        cd_opts <- options$compound_discoverer

        # Filter the consolidated compound table by the user-supplied query.
        all_compounds <- source@compounds
        if (!is.null(cd_opts$compounds_query)) {
            compounds_query <- rlang::parse_expr(cd_opts$compounds_query)
            all_compounds <- all_compounds |> filter(!!compounds_query)
        }

        if (nrow(all_compounds) == 0) {
            stop("No compounds matched the supplied 'compounds_query'.")
        }

        create_compound_chromatograms(
            container, all_compounds, options, cd_opts$rt_extend)
    }
)

#' @rdname create_chromatograms
setMethod(
    f = "create_chromatograms",
    signature = c(source = "LipidSearchSource", features = "NULL"),
    definition = function(source, container, options, features) {
        ls_opts <- options$lipid_search
        all_compounds <- source@compounds

        # The query selects *lipids*, not lipid/sample rows: a row-level
        # predicate such as `grade == "A"` picks the lipids that were graded A
        # somewhere, and every requested sample keeps its extraction row -
        # including samples the result file does not cover.
        if (!is.null(ls_opts$lipids_query)) {
            lipids_query <- rlang::parse_expr(ls_opts$lipids_query)
            matched <- all_compounds |> filter(!!lipids_query)
            all_compounds <- all_compounds |>
                filter(.data$name %in% matched$name)
        }

        if (nrow(all_compounds) == 0) {
            stop("No lipids matched the supplied 'lipids_query'.")
        }

        create_compound_chromatograms(
            container, all_compounds, options, ls_opts$rt_extend)
    }
)

#' @rdname create_chromatograms
setMethod(
    f = "create_chromatograms",
    signature = c(source = "MChromatograms", features = "NULL"),
    definition = function(source, container, options, features) {
        samples <- .select_samples(
            container, options$chromatograms$sample_ids)
        cells <- expand.grid(
            row = seq_len(nrow(source)),
            sample_index = samples$sample_index)

        chrom_data <- tibble(sample_index = cells$sample_index)
        peaks <- .chromatogram_peaks(source, cells)

        container@chromatograms <- as_memory_chromatograms(chrom_data, peaks)
        container@detected_peaks <- tibble()
        container
    }
)

#' @rdname create_chromatograms
setMethod(
    f = "create_chromatograms",
    signature = c(source = "XChromatograms", features = "NULL"),
    definition = function(source, container, options, features) {
        samples <- .select_samples(
            container, options$chromatograms$sample_ids)
        rows <- .xchromatograms_rows(source)
        cells <- expand.grid(
            row = seq_len(nrow(source)),
            sample_index = samples$sample_index)

        chrom_data <- tibble(
            sample_index = cells$sample_index,
            feature_id = rows$feature_id[cells$row],
            feature_mz = rows$feature_mz[cells$row],
            mzMin = rows$mzmin[cells$row],
            mzMax = rows$mzmax[cells$row]
        )
        peaks <- .chromatogram_peaks(source, cells)

        container@chromatograms <- as_memory_chromatograms(chrom_data, peaks)
        container@detected_peaks <- .source_peaks(container, samples)
        container
    }
)

#' @rdname create_chromatograms
setMethod(
    f = "create_chromatograms",
    signature = c(source = "XChromatogram", features = "NULL"),
    definition = function(source, container, options, features) {
        chrom_data <- tibble(sample_index = 1L)
        peaks <- list(data.frame(
            rtime = MSnbase::rtime(source),
            intensity = MSnbase::intensity(source)))

        container@chromatograms <- as_memory_chromatograms(chrom_data, peaks)
        container@detected_peaks <- .source_peaks(container)
        container
    }
)

#' @rdname create_chromatograms
setMethod(
    f = "create_chromatograms",
    signature = c(source = "ANY", features = "NULL"),
    definition = function(source, container, options, features) {
        samples <- .select_samples(
            container, options$chromatograms$sample_ids)
        rt_type <- options$chromatograms$rt_type

        if (rt_type %in% c("corrected", "both") &&
            !has_adjusted_rtime(container)) {
            stop(
                "The data object should be an XcmsExperiment (or XCMSnExp) ",
                "with adjusted retention times to plot the RT adjusted ",
                "chromatograms")
        }

        create_total_chromatograms(
            container,
            samples,
            aggregation_fun = options$chromatograms$aggregation_fun,
            rt_type = rt_type,
            BPPARAM = options$parallel_param)
    }
)

#' @rdname create_chromatograms
setMethod(
    f = "create_chromatograms",
    signature = c(source = "ANY", features = "character"),
    definition = function(source, container, options, features) {
        if (!is_xcms_data(source)) {
            stop(
                "To use feature IDs from the grouped peaks you need to ",
                "provide an xcms object.")
        }

        opts <- options$chromatograms
        samples <- .select_samples(container, opts$sample_ids)
        grouped_peaks <- match_features_by_name(
            get_grouped_peaks(source), features)
        all_peaks <- get_detected_peaks(source)

        # One window per member peak: the peak's m/z with the ppm tolerance and
        # the feature's retention time with the RT tolerance.
        member_peaks <- grouped_peaks |>
            select("name", feature_rt = "rt", "peakidx") |>
            mutate(peak_row = lapply(
                strsplit(.data$peakidx, ","),
                function(x) sort(unique(as.integer(x))))) |>
            select(-"peakidx") |>
            tidyr::unnest("peak_row")

        peaks <- all_peaks[member_peaks$peak_row, , drop = FALSE] |>
            mutate(
                feature_id = member_peaks$name,
                feature_rt = member_peaks$feature_rt) |>
            filter(.data$sample_index %in% samples$sample_index)

        windows <- tibble(
            sample_index = peaks$sample_index,
            feature_id = peaks$feature_id,
            mzMin = peaks$mz - peaks$mz * opts$ppm / 1e6,
            mzMax = peaks$mz + peaks$mz * opts$ppm / 1e6,
            rtMin = peaks$feature_rt - opts$rt_tol,
            rtMax = peaks$feature_rt + opts$rt_tol
        )

        container <- extract_chromatograms(
            container,
            windows,
            adjusted = has_adjusted_rtime(container),
            fill_gaps = opts$fill_gaps,
            BPPARAM = options$parallel_param)
        container@detected_peaks <- normalise_peaks(
            peaks |> select(-"feature_rt"))
        container
    }
)

#' @rdname create_chromatograms
setMethod(
    f = "create_chromatograms",
    signature = c(source = "ANY", features = "ANY"),
    definition = function(source, container, options, features) {
        opts <- options$chromatograms
        samples <- .select_samples(container, opts$sample_ids)
        adjusted <- has_adjusted_rtime(container)

        # Windows without an RT range span the sample's whole acquisition,
        # which is only known from the raw data.
        rt_ranges <- NULL
        if (.needs_full_rt_range(features)) {
            container <- ensure_spectra(container, options$parallel_param)
            sps <- select_rtime(raw_spectra(container), adjusted)
            rt_ranges <- lapply(
                split(Spectra::rtime(sps), sps$sample_index), range)
        }

        windows <- bind_rows(lapply(seq_len(nrow(samples)), function(i) {
            sample_metadata <- samples[i, ]
            feats <- get_features(
                options,
                sample_metadata,
                full_rt_range = rt_ranges[[
                    as.character(sample_metadata$sample_index)]])

            tibble(
                sample_index = rep(
                    sample_metadata$sample_index, length(feats)),
                feature_id = vapply(feats, `[[`, character(1), "feature_id"),
                mzMin = vapply(feats, function(f) f$mzr[1], numeric(1)),
                mzMax = vapply(feats, function(f) f$mzr[2], numeric(1)),
                rtMin = vapply(feats, function(f) f$rtr[1], numeric(1)),
                rtMax = vapply(feats, function(f) f$rtr[2], numeric(1))
            )
        }))

        container <- extract_chromatograms(
            container,
            windows,
            adjusted = adjusted,
            fill_gaps = opts$fill_gaps,
            BPPARAM = options$parallel_param)
        container@detected_peaks <- peaks_in_windows(
            get_detected_peaks(source), windows)
        container
    }
)

# Feature specifications that give no RT range extract over the whole run.
.needs_full_rt_range <- function(features) {
    features <- as.data.frame(features)
    has_rt <- "rt" %in% colnames(features) ||
        all(c("rtmin", "rtmax") %in% colnames(features))
    if (!has_rt) {
        return(TRUE)
    }
    rt_cols <- intersect(c("rt", "rtmin", "rtmax"), colnames(features))
    anyNA(features[, rt_cols])
}

#' Create base peak, total ion, or mean ion chromatograms
#'
#' One chromatogram per sample over all its MS1 scans. The scan header values
#' are used when the format reports them reliably, so no peaks are read;
#' netCDF files and in-memory data without header values are summarised from
#' their peaks.
#'
#' @param container An `lcmsPlotDataContainer`.
#' @param samples A `data.frame` of the samples to process.
#' @param aggregation_fun One of `"max"` (base peak), `"sum"` (total ion
#' current), or `"mean"` (mean intensity of the scan's peaks).
#' @param rt_type One of `"uncorrected"`, `"corrected"`, or `"both"`.
#' @param BPPARAM A `BiocParallelParam`, or `NULL` to process sequentially.
#' @return The `lcmsPlotDataContainer` with its `chromatograms` slot set and
#' no detected peaks.
#' @keywords internal
create_total_chromatograms <- function(
    container,
    samples,
    aggregation_fun,
    rt_type = "uncorrected",
    BPPARAM = NULL
) {
    aggregation_fun <- match.arg(aggregation_fun, c("max", "sum", "mean"))
    paths <- MsExperiment::sampleData(container@experiment)$sample_path
    use_rawrr <- .use_rawrr(paths)

    if (!all(use_rawrr[samples$sample_index])) {
        container <- ensure_spectra(container, BPPARAM)
    }

    rt_types <- switch(
        rt_type,
        uncorrected = NA_character_,
        corrected = "RT corrected",
        both = c("RT uncorrected", "RT corrected"))

    results <- .apply_samples(samples$sample_index, function(s) {
        if (use_rawrr[s]) {
            peaks <- rawrr_bpc_tic(paths[s], aggregation_fun)
            return(list(
                chrom_data = tibble(sample_index = s, rt_type = rt_types),
                peaks = list(peaks)))
        }

        sps <- raw_spectra(container, s, ms_level = 1L)
        intensity <- .scan_intensity(sps, aggregation_fun)
        peaks <- lapply(rt_types, function(type) {
            rtime <- if (identical(type, "RT corrected")) {
                sps$rtime_adjusted
            } else {
                Spectra::rtime(sps)
            }
            data.frame(rtime = rtime, intensity = intensity)
        })

        list(
            chrom_data = tibble(sample_index = s, rt_type = rt_types),
            peaks = peaks)
    }, BPPARAM)

    container@chromatograms <- as_memory_chromatograms(
        bind_rows(lapply(results, `[[`, "chrom_data")),
        unlist(lapply(results, `[[`, "peaks"), recursive = FALSE))
    container@detected_peaks <- tibble()
    container
}

#' Create chromatograms for compound tables
#'
#' Shared by the compound-centric data sources (Compound Discoverer scripting
#' node, LipidSearch). Each row of `all_compounds` describes where to extract a
#' compound in one sample: an `mz`, an `rtmin`/`rtmax` window, and the
#' `sample_index` it belongs to.
#'
#' Compounds with no `mz` or no RT window are skipped, as are compounds absent
#' from a given sample. When `all_compounds` carries a `detected` column, a
#' detected peak - the one the `highlight_peaks` layer draws - is kept only
#' where it is `TRUE`, so compounds that were looked for but not found are
#' still plotted without a degenerate highlight.
#'
#' @param container An `lcmsPlotDataContainer`.
#' @param all_compounds A `data.frame` with at least `name`, `mz`, `rt`,
#' `rtmin`, `rtmax`, `into`, `maxo`, and `sample_index` columns. `formula`,
#' `adduct`, `class`, `sub_class`, and `detected` are used when present.
#' @param options A `list` of plot options.
#' @param rt_extend A `numeric` value indicating how much (in seconds) to extend
#' the RT window on each side of the compound peak.
#' @return The `lcmsPlotDataContainer` with its `chromatograms` and
#' `detected_peaks` slots set.
#' @keywords internal
create_compound_chromatograms <- function(
    container,
    all_compounds,
    options,
    rt_extend
) {
    opts <- options$chromatograms
    samples <- .select_samples(container, opts$sample_ids)
    annotation_cols <- intersect(
        c("name", "formula", "adduct", "mz", "class", "sub_class"),
        colnames(all_compounds))
    compound_names <- unique(all_compounds$name)

    entries <- all_compounds |>
        filter(.data$sample_index %in% samples$sample_index) |>
        distinct(.data$sample_index, .data$name, .keep_all = TRUE) |>
        filter(!is.na(.data$mz), !is.na(.data$rtmin), !is.na(.data$rtmax)) |>
        arrange(
            match(.data$sample_index, samples$sample_index),
            match(.data$name, compound_names))

    windows <- entries |>
        select(all_of(c("sample_index", annotation_cols))) |>
        mutate(
            feature_id = .data$name,
            mzMin = entries$mz - entries$mz * opts$ppm / 1e6,
            mzMax = entries$mz + entries$mz * opts$ppm / 1e6,
            rtMin = entries$rtmin - rt_extend,
            rtMax = entries$rtmax + rt_extend
        )

    container <- extract_chromatograms(
        container,
        windows,
        fill_gaps = opts$fill_gaps,
        BPPARAM = options$parallel_param)

    # Sources that look compounds up in samples where they were never reported
    # flag the difference with a `detected` column; those entries have no peak
    # to highlight.
    detected <- if ("detected" %in% colnames(entries)) {
        entries |> filter(.data$detected %in% TRUE)
    } else {
        entries
    }

    container@detected_peaks <- normalise_peaks(
        detected |>
            mutate(feature_id = .data$name) |>
            select(all_of(c(
                "name", "feature_id", "mz", "rt", "rtmin", "rtmax", "into",
                "maxo", "sample_index"))))
    container
}

# The retention times and intensities of the chromatograms held in an
# MChromatograms/XChromatograms object, for the given (row, sample) cells.
.chromatogram_peaks <- function(source, cells) {
    lapply(seq_len(nrow(cells)), function(k) {
        chrom <- source[cells$row[k], cells$sample_index[k]]
        data.frame(
            rtime = MSnbase::rtime(chrom),
            intensity = MSnbase::intensity(chrom))
    })
}

# One feature identifier, m/z window, and representative m/z per EIC row.
.xchromatograms_rows <- function(source) {
    n_rows <- nrow(source)

    # mz(<XChromatograms>) is an mzmin/mzmax matrix, one row per EIC, which is
    # what orders the series when lp_chromatogram(stacked = ) is used.
    mz_windows <- MSnbase::mz(source)

    # featureChromatograms() records which row each feature came from.
    feat_defs <- xcms::featureDefinitions(source)
    feat_rows <- if ("row" %in% colnames(feat_defs)) feat_defs$row else integer()
    row_names <- get_feature_names(source)

    feature_id <- vapply(seq_len(n_rows), function(i) {
        rts <- MSnbase::rtime(source[i, 1L])
        mz_center <- mean(mz_windows[i, ])
        rt_center <- if (length(rts) > 0L) mean(range(rts)) else NA_real_
        row_feature_names <- rownames(feat_defs)[feat_rows == i]

        if (!is.null(row_names) && !is.na(row_names[i])) {
            row_names[i]
        } else if (length(row_feature_names) == 1L) {
            row_feature_names
        } else if (!is.na(rt_center)) {
            sprintf("M%dT%d", round(mz_center), round(rt_center))
        } else {
            sprintf("M%d", round(mz_center))
        }
    }, character(1))

    tibble(
        feature_id = feature_id,
        feature_mz = rowMeans(matrix(as.numeric(mz_windows), nrow = n_rows)),
        mzmin = as.numeric(mz_windows[, 1L]),
        mzmax = as.numeric(mz_windows[, 2L])
    )
}

# Drop the data points without an intensity from every chromatogram.
drop_missing_intensities <- function(chromatograms) {
    if (length(chromatograms) == 0) {
        return(chromatograms)
    }
    peaks <- lapply(
        Chromatograms::peaksData(chromatograms),
        function(p) p[!is.na(p$intensity), , drop = FALSE])
    as_memory_chromatograms(
        as.data.frame(Chromatograms::chromData(chromatograms)), peaks)
}
