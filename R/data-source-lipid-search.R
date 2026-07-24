# Lipid-level (sample-independent) columns of a LipidSearch result file.
# Each logical field is resolved against a list of candidate
# column names so that differences between LipidSearch versions can
# be absorbed by adding a candidate. Columns absent from a given version resolve
# to NA.
.ls_lipid_cols <- list(
    rej = c("Rej.", "Rej"),
    lipid_ion = c("LipidIon", "LipidID"),
    lipid_group = c("LipidGroup"),
    class = c("Class"),
    sub_class = c("SubClass"),        # 5.2 only
    fatty_acid = c("FattyAcid"),      # 4.2 only
    calc_mz = c("CalcMz"),
    formula = c("IonFormula"),
    adduct = c("AdductIon"),          # 5.2 only (explicit adduct column)
    base_rt = c("BaseRt")             # 5.2 only (per-lipid consensus RT)
)

# Lipid-level fields that are numeric rather than character. `rej` is handled
# separately (it is 0/1 in 4.2 but false/true in 5.2).
.ls_lipid_numeric <- c("calc_mz", "base_rt")

# Per-sample column base names. In the file each of these is suffixed with the
# per-injection key in square brackets, e.g. "Area[c-1]" (4.2), "Area[s1-1]"
# (5.2).
.ls_sample_cols <- list(
    into = c("Area"),
    maxo = c("Height"),
    rt = c("Rt", "ObsRt"),
    obs_mz = c("ObsMz"),
    hwhm_l = c("Hwhm(L)"),
    hwhm_r = c("Hwhm(R)"),
    grade = c("Grade"),
    mscore = c("mScore", "IDScore"),
    sn = c("S/N")
)

# Per-sample fields that are numeric rather than character.
.ls_sample_numeric <- c(
    "into", "maxo", "rt", "obs_mz", "hwhm_l", "hwhm_r", "mscore", "sn"
)

# Sample declaration lines (4.2) look like: #[c-1]:PT_C18P_..._QC.raw
.ls_sample_line <- "^#\\[([^]]+)\\]:(.*)$"

# 5.2 group-level column bases whose bracket keys define the sample keys
# (s1, s2, ...). The first present is used.
.ls_group_cols <- c("OrgMeanArea", "MeanArea", "OrgMeanHeight")

#' Extract the bracket key from a suffixed column name
#'
#' `"Area[s1-1]"` -> `"s1-1"`; a name without a bracket suffix yields `NA`.
#'
#' @param cols A `character` vector of column names.
#' @return A `character` vector of bracket keys.
#' @keywords internal
.ls_bracket_key <- function(cols) {
    ifelse(
        grepl("\\[[^]]+\\]$", cols),
        sub("^.*\\[([^]]+)\\]$", "\\1", cols),
        NA_character_
    )
}

#' Discover the 5.2 sample declarations from the table columns
#'
#' LipidSearch 5.2 files carry no `#[key]:file` declarations. The sample keys are
#' taken from the group-level `OrgMeanArea[...]` columns (`s1`, `s2`, ...); each
#' key's per-injection data lives in columns suffixed `<key>-<n>` (e.g.
#' `Area[s1-1]`), so the first such injection key is resolved per sample.
#'
#' @param cols A `character` vector of the table's column names.
#' @return A `tibble` of `sample_key`, `sample_group`, `injection_key` and
#' `file_name` (always `NA`), or `NULL` when no group columns are found.
#' @keywords internal
.ls_discover_samples <- function(cols) {
    group_base <- .ls_group_cols[
        vapply(
            .ls_group_cols,
            function(b) any(grepl(paste0("^", b, "\\[[^]]+\\]$"), cols)),
            logical(1)
        )
    ]
    if (length(group_base) == 0) return(NULL)

    group_cols <- grep(
        paste0("^", group_base[[1]], "\\[[^]]+\\]$"), cols, value = TRUE)
    keys <- .ls_bracket_key(group_cols)

    all_keys <- .ls_bracket_key(cols)
    injection_key <- vapply(keys, function(k) {
        hit <- all_keys[
            !is.na(all_keys) & grepl(paste0("^", .ls_escape(k), "-[0-9]+$"),
                                     all_keys)
        ]
        if (length(hit) > 0) hit[[1]] else k
    }, character(1))

    tibble(
        sample_key = keys,
        sample_group = keys,
        injection_key = injection_key,
        file_name = NA_character_
    )
}

#' Escape a string for literal use inside a regular expression
#' @keywords internal
.ls_escape <- function(x) {
    gsub("([].|()\\^{}+$*?[])", "\\\\\\1", x)
}

#' Read a LipidSearch result file
#'
#' Reads a LipidSearch result file and detects its version. LipidSearch 4.2 files
#' begin with a block of `#[key]:file` sample declarations and `#key:value`
#' settings, followed by the tab-delimited lipid table. LipidSearch 5.2 files
#' carry no declarations; the table starts on the first line and the samples are
#' discovered from the `OrgMeanArea[...]` columns.
#'
#' @param path A `character` path to a LipidSearch result file.
#' @return A named `list` with elements `version` (`"4.2"` or `"5.2"`), `samples`
#' (a `tibble` of `sample_key`, `sample_group`, `injection_key` and `file_name`),
#' `settings` (a named `list`), and `data` (a `tibble` of the lipid table).
#' @keywords internal
.ls_read_results <- function(path) {
    if (!is.character(path) || length(path) != 1 || !file.exists(path)) {
        stop("results_path must be a path to an existing LipidSearch file.")
    }

    lines <- readLines(path, warn = FALSE)

    comment_lines <- lines[startsWith(lines, "#")]
    decl_lines <- grep(.ls_sample_line, comment_lines, value = TRUE)
    is_v4 <- length(decl_lines) > 0

    # The lipid table starts at the first non-comment, non-blank line;
    # everything from there on is a tab-delimited table with a header.
    is_data <- !startsWith(lines, "#") & nzchar(trimws(lines))
    if (!any(is_data)) {
        stop("No lipid table found in '", path, "'.")
    }

    data <- utils::read.table(
        text = lines[seq(which(is_data)[[1]], length(lines))],
        sep = "\t",
        header = TRUE,
        check.names = FALSE,
        stringsAsFactors = FALSE,
        quote = "",
        comment.char = ""
    )

    # The 4.2 header row ends with a trailing tab, yielding one unnamed column
    # that has to go before the table can become a tibble.
    data <- as_tibble(data[, nzchar(colnames(data)), drop = FALSE])

    if (is_v4) {
        keys <- sub(.ls_sample_line, "\\1", decl_lines)
        samples <- tibble(
            sample_key = keys,
            # "c-1" -> group "c". Keys without a "-" are their own group.
            sample_group = sub("-[^-]*$", "", keys),
            # In 4.2 the per-injection columns are keyed by the sample key
            # itself (e.g. Area[c-1]).
            injection_key = keys,
            file_name = trimws(sub(.ls_sample_line, "\\2", decl_lines))
        )

        setting_lines <- setdiff(comment_lines, decl_lines)
        setting_lines <- setting_lines[grepl(":", setting_lines, fixed = TRUE)]
        settings <- as.list(trimws(sub("^[^:]*:", "", setting_lines)))
        names(settings) <- trimws(sub("^#([^:]*):.*$", "\\1", setting_lines))

        return(list(
            version = "4.2", samples = samples,
            settings = settings, data = data))
    }

    samples <- .ls_discover_samples(colnames(data))
    if (is.null(samples)) {
        stop(
            "Could not identify any samples in '", path, "'. Expected either ",
            "'#[c-1]:sample.raw' declarations (LipidSearch 4.2) or per-sample ",
            "'OrgMeanArea[...]' columns (LipidSearch 5.2)."
        )
    }

    list(version = "5.2", samples = samples, settings = list(), data = data)
}

#' Coerce a LipidSearch column to the expected type
#'
#' Cells for undetected lipids are left empty, which makes `read.table()` type
#' otherwise-numeric columns as `character`.
#'
#' @param values The raw column values, or `NULL` when the column is absent.
#' @param numeric A `logical` value indicating whether the field is numeric.
#' @param n The number of rows to return when `values` is `NULL`.
#' @return A `numeric` or `character` vector of length `n`.
#' @keywords internal
.ls_coerce <- function(values, numeric, n) {
    if (is.null(values)) {
        return(rep(if (numeric) NA_real_ else NA_character_, n))
    }

    if (numeric) {
        suppressWarnings(as.numeric(values))
    } else {
        values <- as.character(values)
        ifelse(nzchar(trimws(values)), values, NA_character_)
    }
}

#' Resolve the lipid-level columns of a LipidSearch table
#'
#' @param df The lipid table from `.ls_read_results()`.
#' @return A `tibble` with the standardised lipid-level columns. Fields absent
#' from the export are filled with `NA`.
#' @keywords internal
.ls_lipid_fields <- function(df) {
    resolved <- vapply(
        .ls_lipid_cols, function(x) first_matching_column(df, x), character(1))

    required <- c("lipid_ion", "calc_mz")
    missing <- required[is.na(resolved[required])]
    if (length(missing) > 0) {
        stop(
            "Missing required LipidSearch columns: ",
            paste(
                vapply(.ls_lipid_cols[missing], `[[`, character(1), 1),
                collapse = ", "
            )
        )
    }

    fields <- lapply(names(.ls_lipid_cols), function(field) {
        col <- resolved[[field]]
        values <- if (is.na(col)) NULL else df[[col]]

        if (field == "rej") {
            .ls_parse_rejected(values, nrow(df))
        } else {
            .ls_coerce(values, field %in% .ls_lipid_numeric, nrow(df))
        }
    })

    names(fields) <- names(.ls_lipid_cols)
    as_tibble(fields)
}

#' Parse the LipidSearch rejection flag as a logical
#'
#' The flag is `0`/`1` in LipidSearch 4.2 and `false`/`true` in 5.2. Absent or
#' unparseable values are treated as *not* rejected.
#'
#' @param values The raw column values, or `NULL` when the column is absent.
#' @param n The number of rows to return when `values` is `NULL`.
#' @return A `logical` vector of length `n`, never `NA`.
#' @keywords internal
.ls_parse_rejected <- function(values, n) {
    if (is.null(values)) {
        return(rep(FALSE, n))
    }
    v <- tolower(trimws(as.character(values)))
    v %in% c("1", "true", "t", "yes", "y")
}

#' Extract one sample's columns from a LipidSearch table
#'
#' @param df The lipid table from `.ls_read_results()`.
#' @param injection_key A `character` value giving the per-injection column key,
#' e.g. `"c-1"` (4.2) or `"s1-1"` (5.2).
#' @return A `tibble` with one row per lipid and the standardised per-sample
#' columns. Fields absent from the export are filled with `NA`. In 5.2 the
#' undetected marker `"NP"` becomes `NA` for numeric fields.
#' @keywords internal
.ls_sample_fields <- function(df, injection_key) {
    fields <- lapply(names(.ls_sample_cols), function(field) {
        col <- first_matching_column(
            df,
            paste0(.ls_sample_cols[[field]], "[", injection_key, "]")
        )
        .ls_coerce(
            if (is.na(col)) NULL else df[[col]],
            field %in% .ls_sample_numeric,
            nrow(df)
        )
    })

    names(fields) <- names(.ls_sample_cols)
    as_tibble(fields)
}

#' Reshape a LipidSearch table into one row per lipid and declared sample
#'
#' @param res The parsed result of `.ls_read_results()`.
#' @param keep_rejected A `logical` value indicating whether rows flagged as
#' rejected should be retained.
#' @return A `tibble` with one row per (lipid row, declared sample), carrying
#' the lipid-level fields plus `lipid_row` (the source row index), `sample_key`
#' and `detected`.
#' @keywords internal
.ls_build_lipids <- function(res, keep_rejected) {
    lipids <- .ls_lipid_fields(res$data) |>
        mutate(lipid_row = row_number())

    if (!keep_rejected) {
        lipids <- lipids |> filter(!.data$rej)
    }

    if (nrow(lipids) == 0) {
        stop("No lipids left after filtering the LipidSearch result file.")
    }

    rows <- lapply(seq_len(nrow(res$samples)), function(i) {
        sample <- res$samples[i, ]
        per_sample <- .ls_sample_fields(
            res$data, sample$injection_key)[lipids$lipid_row, ]
        bind_cols(lipids, per_sample) |>
            mutate(
                sample_key = sample$sample_key,
                # A missing/zero area means LipidSearch reported the lipid on
                # this row but did not integrate a peak in this sample.
                detected = !is.na(.data$into) & .data$into > 0
            )
    })

    bind_rows(rows)
}

#' Synthesize unique lipid labels
#'
#' The plot keys compound identity on `name`, but `LipidIon` is not unique: the
#' same ion is commonly reported at several retention times. Duplicated ions are
#' disambiguated with their consensus RT, falling back to the source row index.
#'
#' @param lipid_ion A `character` vector of lipid ion names.
#' @param rt A `numeric` vector of consensus retention times, in minutes.
#' @param lipid_row An `integer` vector of source row indices.
#' @return A `character` vector of unique labels.
#' @keywords internal
.ls_lipid_labels <- function(lipid_ion, rt, lipid_row) {
    base <- ifelse(
        !is.na(lipid_ion) & nzchar(trimws(lipid_ion)),
        lipid_ion,
        paste0("Lipid ", lipid_row)
    )

    dup <- base %in% base[duplicated(base)]
    labelled <- ifelse(
        dup & !is.na(rt),
        sprintf("%s [RT %.2f]", base, rt),
        base
    )

    still_dup <- labelled %in% labelled[duplicated(labelled)]
    ifelse(still_dup, paste0(labelled, " [#", lipid_row, "]"), labelled)
}

#' Derive the adduct from a LipidSearch ion name
#'
#' LipidSearch ion names carry the adduct as a suffix, e.g. `AEA(18:2)+H` or
#' `PC(16:0_18:1)+NH4`.
#'
#' @param lipid_ion A `character` vector of lipid ion names.
#' @return A `character` vector of adducts, `NA` where none can be determined.
#' @keywords internal
.ls_adduct <- function(lipid_ion) {
    adduct <- sub("^.*\\)", "", lipid_ion)
    ifelse(
        !is.na(adduct) & nzchar(adduct) & adduct != lipid_ion,
        adduct,
        NA_character_
    )
}

#' Compute per-lipid consensus extraction values
#'
#' These drive chromatogram extraction in samples that are absent from the
#' result file, and act as a fallback for declared samples in which the lipid
#' was not detected. Consensus values are taken over the samples where the
#' lipid *was* detected, falling back to all declared samples, the reported
#' `BaseRt` (5.2), and finally the theoretical m/z.
#'
#' @param lipids A `tibble` as returned by `.ls_build_lipids()`.
#' @return A `tibble` with one row per `lipid_row` and columns `lipid_row`,
#' `mz`, `rt`, `rtmin` and `rtmax` (retention times in minutes).
#' @keywords internal
.ls_consensus <- function(lipids) {
    # median() of an all-NA vector with na.rm = TRUE is NA, which is exactly
    # the fallback signal the coalesce() chains below rely on.
    med <- function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)

    lipids |>
        group_by(.data$lipid_row) |>
        summarise(
            calc_mz = dplyr::first(.data$calc_mz),
            base_rt = dplyr::first(.data$base_rt),
            mz_detected = med(.data$obs_mz[.data$detected]),
            mz_any = med(.data$obs_mz),
            rt_detected = med(.data$rt[.data$detected]),
            rt_any = med(.data$rt),
            hwhm_l = med(.data$hwhm_l[.data$detected]),
            hwhm_r = med(.data$hwhm_r[.data$detected]),
            .groups = "drop"
        ) |>
        mutate(
            mz = dplyr::coalesce(
                .data$mz_detected, .data$mz_any, .data$calc_mz),
            # Prefer detected RTs; otherwise 5.2's per-lipid BaseRt lets even a
            # lipid detected in no sample still be extracted.
            rt = dplyr::coalesce(
                .data$rt_detected, .data$base_rt, .data$rt_any),
            rtmin = .data$rt - dplyr::coalesce(.data$hwhm_l, 0),
            rtmax = .data$rt + dplyr::coalesce(.data$hwhm_r, 0)
        ) |>
        select(
            .data$lipid_row, .data$mz, .data$rt, .data$rtmin, .data$rtmax)
}

#' Match sample paths to the file's declared samples
#'
#' Returns, for each `sample_paths` entry, the index into `declared` it maps to
#' (`NA` for an extra sample). LipidSearch 4.2 files name their raw files, so
#' paths are matched by basename without extension. LipidSearch 5.2 files carry
#' no filenames: a named `sample_paths` vector is matched by its names against
#' the sample keys (`s1`, `s2`, ...), and an unnamed vector is matched
#' positionally.
#'
#' @param declared A `tibble` of sample declarations from `.ls_read_results()`.
#' @param sample_paths A `character` vector of raw sample file paths.
#' @return An `integer` vector of `nrow == length(sample_paths)` of row indices
#' into `declared`, or `NA` for extra samples.
#' @keywords internal
.ls_match_samples <- function(declared, sample_paths) {
    has_files <- any(!is.na(declared$file_name))

    if (has_files) {
        # 4.2: match by basename without extension.
        key <- function(x) tools::file_path_sans_ext(basename(x))
        return(match(key(sample_paths), key(declared$file_name)))
    }

    # 5.2: match named entries by sample key, otherwise positionally.
    nm <- names(sample_paths)
    if (!is.null(nm) && any(nzchar(nm))) {
        idx <- match(nm, declared$sample_key)
        idx[!nzchar(nm)] <- NA_integer_
        return(idx)
    }

    ifelse(
        seq_along(sample_paths) <= nrow(declared),
        seq_along(sample_paths),
        NA_integer_
    )
}

#' Build the sample metadata for a LipidSearch source
#'
#' `sample_paths` is authoritative: every supplied path becomes a sample. Paths
#' are matched to the result file's samples by `.ls_match_samples()` - by
#' filename for 4.2, by sample key / position for 5.2. Paths that match nothing
#' are kept as extra samples, plotted from each lipid's consensus m/z and RT
#' window.
#'
#' @param declared A `tibble` of sample declarations from `.ls_read_results()`.
#' @param sample_paths A `character` vector of raw sample file paths.
#' @param metadata A `data.frame` (optional) of additional sample metadata with
#' one row per `sample_paths` entry.
#' @return A `tibble` with `sample_index`, `sample_id`, `sample_path`,
#' `sample_key`, `sample_group` and `in_results`.
#' @keywords internal
.ls_build_metadata <- function(declared, sample_paths, metadata = NULL) {
    if (!is.character(sample_paths) || length(sample_paths) == 0) {
        stop("sample_paths must be a non-empty character vector.")
    }

    idx <- .ls_match_samples(declared, sample_paths)

    unmatched <- setdiff(seq_len(nrow(declared)), idx)
    if (length(unmatched) > 0) {
        labels <- ifelse(
            is.na(declared$file_name[unmatched]),
            declared$sample_key[unmatched],
            declared$file_name[unmatched]
        )
        message(
            "LipidSearch: no sample path supplied for ",
            length(unmatched), " declared sample(s): ",
            paste(labels, collapse = ", ")
        )
    }

    sample_metadata <- tibble(
        sample_index = seq_along(sample_paths),
        sample_id = tools::file_path_sans_ext(basename(sample_paths)),
        sample_path = unname(sample_paths),
        sample_key = declared$sample_key[idx],
        sample_group = declared$sample_group[idx],
        in_results = !is.na(idx)
    )

    if (!is.null(metadata)) {
        if (nrow(metadata) != length(sample_paths)) {
            stop("metadata must have one row per sample path.")
        }
        sample_metadata <- as_tibble(
            cbind(sample_metadata, as_tibble(metadata))
        )
    }

    sample_metadata
}

#' LipidSearch data source
#'
#' An S4 class wrapping the lipid annotations exported by LipidSearch.
#'
#' @slot name A `character` value identifying the data source.
#' @slot version A `character` value giving the detected LipidSearch format,
#' `"4.2"` or `"5.2"`.
#' @slot metadata A `data.frame` of sample metadata (one row per supplied
#' sample path), including a `sample_path` column with the raw file paths and an
#' `in_results` flag marking samples declared in the result file.
#' @slot peaks A `data.frame` of detected peaks (one row per lipid and sample in
#' which LipidSearch integrated a peak) with the standard `mz`, `rt`, `rtmin`,
#' `rtmax`, `into`, `maxo` and `sample_index` columns, plus `name`, `formula`,
#' `adduct`, `class` and `grade` annotations.
#' @slot compounds A `data.frame` with one row per lipid and sample - including
#' samples absent from the result file - used to drive chromatogram extraction.
#' @export
setClass(
    "LipidSearchSource",
    slots = list(
        name = "character",
        version = "character",
        metadata = "data.frame",
        peaks = "data.frame",
        compounds = "data.frame"
    )
)

#' Show a summary of a `LipidSearchSource` object
#'
#' @param object An instance of class `LipidSearchSource`.
#' @return Invisible \code{NULL}
#' @export
setMethod(
    f = "show",
    signature = "LipidSearchSource",
    function(object) {
        cat("Object of class", class(object), "\n")
        cat(" Name:", object@name, "\n")
        cat(" Version:", object@version, "\n")
        cat(
            " Samples:", nrow(object@metadata),
            paste0("(", sum(object@metadata$in_results), " in results)"), "\n",
            " Lipids:", length(unique(object@compounds$name)), "\n",
            " Peaks:", nrow(object@peaks), "rows\n"
        )
    }
)

#' Create a data source from a LipidSearch result file
#'
#' Reads a LipidSearch result file and builds a data source that can be plotted
#' with `lcmsPlot()`. Both **LipidSearch 4.2** and **5.2** exports are supported;
#' the version is auto-detected. Chromatograms are later extracted from the raw
#' files supplied via `sample_paths`.
#'
#' @section Sample-to-file mapping:
#' A LipidSearch 4.2 file names the analysed raw files in its header
#' (`#[c-1]:sample.raw`), so paths are matched to samples by file basename
#' without extension - a result file that lists ThermoFisher `.raw` files works
#' with converted `.mzML` files. A LipidSearch 5.2 file carries **no** filenames;
#' its samples are keyed `s1`, `s2`, ... (from the `OrgMeanArea[...]` columns).
#' For 5.2, supply either a **named** `sample_paths` vector whose names are those
#' keys (`c("s1" = "a.mzML", "s2" = "b.mzML")`, order-independent, any subset),
#' or an **unnamed** vector matched positionally (`s1`, `s2`, ... in order).
#'
#' @section Samples beyond the result file:
#' `sample_paths` determines which samples are plotted.
#' A path that matches no sample (an unnamed extra for 5.2, or an unmatched
#' basename for 4.2) is still a full sample: for every queried lipid its
#' chromatogram is extracted using the lipid's consensus m/z and retention-time
#' window (the median across the samples in which LipidSearch did detect it,
#' falling back to the reported `BaseRt` for 5.2), it simply has no reported
#' peak to highlight.
#'
#' @param results_path A `character` path to a LipidSearch result file (4.2 or
#' 5.2).
#' @param sample_paths A `character` vector of paths to the raw sample files
#' (e.g. `.raw` or `.mzML`) used to extract chromatograms. For 5.2 files this may
#' be named by sample key (see *Sample-to-file mapping*).
#' @param metadata A `data.frame` (optional) of additional sample metadata. It
#' must contain one row per sample, aligned with `sample_paths`. When supplied,
#' its columns are column-bound onto the generated metadata.
#' @param keep_rejected A `logical` value indicating whether lipids flagged as
#' rejected (`Rej.`/`Rej`) should be retained. Defaults to `FALSE`.
#' @return An object of class `LipidSearchSource`.
#' @seealso [lp_lipid_search()]
#' @export
#' @examples
#' \dontrun{
#' ## LipidSearch 4.2 (raw files declared in the header, matched by basename):
#' ds <- LipidSearchSource(
#'     results_path = "LIPIDS_POS_5ppm.txt",
#'     sample_paths = c("path/to/sample1.raw", "path/to/extra_sample.raw")
#' )
#'
#' ## LipidSearch 5.2 (no filenames; map by sample key, extra sample unnamed):
#' ds <- LipidSearchSource(
#'     results_path = "LIPIDS_POS_5ppm_LS5.txt",
#'     sample_paths = c(
#'         "s1" = "path/to/sample1.mzML",
#'         "s2" = "path/to/sample2.mzML",
#'         "path/to/extra_sample.mzML"
#'     )
#' )
#'
#' lcmsPlot(ds) +
#'     lp_lipid_search(
#'         lipids_query = 'class == "AcCa" & grade == "A"',
#'         rt_extend = 30
#'     ) +
#'     lp_chromatogram(highlight_peaks = TRUE) +
#'     lp_grid(rows = "sample_id", cols = "name", free_x = TRUE)
#' }
LipidSearchSource <- function(
    results_path,
    sample_paths,
    metadata = NULL,
    keep_rejected = FALSE
) {
    res <- .ls_read_results(results_path)
    lipids <- .ls_build_lipids(res, keep_rejected)
    sample_metadata <- .ls_build_metadata(res$samples, sample_paths, metadata)

    # Sample-independent annotations, plus the consensus extraction window used
    # for samples that the result file does not cover.
    annotations <- lipids |>
        distinct(.data$lipid_row, .keep_all = TRUE) |>
        select(
            .data$lipid_row, .data$lipid_ion, .data$lipid_group, .data$class,
            .data$sub_class, .data$fatty_acid, .data$formula, .data$calc_mz,
            .data$rej, adduct_col = .data$adduct
        ) |>
        left_join(.ls_consensus(lipids), by = "lipid_row") |>
        mutate(
            name = .ls_lipid_labels(
                .data$lipid_ion, .data$rt, .data$lipid_row),
            # The explicit AdductIon column (5.2) wins; otherwise derive the
            # adduct from the ion-name suffix (4.2).
            adduct = dplyr::coalesce(
                .data$adduct_col, .ls_adduct(.data$lipid_ion))
        ) |>
        select(-.data$adduct_col) |>
        rename(
            consensus_mz = "mz",
            consensus_rt = "rt",
            consensus_rtmin = "rtmin",
            consensus_rtmax = "rtmax"
        )

    measurements <- lipids |>
        select(
            .data$lipid_row, .data$sample_key, .data$into, .data$maxo,
            .data$obs_mz, sample_rt = .data$rt, .data$hwhm_l, .data$hwhm_r,
            .data$grade, .data$mscore, .data$sn, .data$detected
        )

    # One extraction row per lipid and per *requested* sample. Samples the
    # result file declares contribute their own values; the rest - and declared
    # samples in which nothing was found - fall back to the consensus.
    compounds <- sample_metadata |>
        select(.data$sample_index, .data$sample_key) |>
        cross_join(annotations) |>
        left_join(measurements, by = c("lipid_row", "sample_key")) |>
        mutate(
            detected = !is.na(.data$detected) & .data$detected,
            mz = dplyr::coalesce(.data$obs_mz, .data$consensus_mz),
            rt = dplyr::coalesce(.data$sample_rt, .data$consensus_rt),
            rtmin = ifelse(
                is.na(.data$sample_rt),
                .data$consensus_rtmin,
                .data$sample_rt - dplyr::coalesce(.data$hwhm_l, 0)
            ),
            rtmax = ifelse(
                is.na(.data$sample_rt),
                .data$consensus_rtmax,
                .data$sample_rt + dplyr::coalesce(.data$hwhm_r, 0)
            ),
            # Only integrated peaks carry an area/height; the zeros reported for
            # undetected lipids would otherwise draw degenerate highlights.
            into = ifelse(.data$detected, .data$into, NA_real_),
            maxo = ifelse(.data$detected, .data$maxo, NA_real_)
        ) |>
        select(
            .data$name, .data$lipid_ion, .data$lipid_group, .data$class,
            .data$sub_class, .data$fatty_acid, .data$formula, .data$adduct,
            .data$calc_mz, .data$rej, .data$mz, .data$rt, .data$rtmin,
            .data$rtmax, .data$into, .data$maxo, .data$grade, .data$mscore,
            .data$sn, .data$detected, .data$sample_index, .data$lipid_row
        ) |>
        # LipidSearch reports retention times in minutes.
        convert_rt_to_seconds() |>
        add_compound_rank(rank_column = "lipid_rank")

    peaks <- compounds |>
        filter(.data$detected) |>
        select(
            .data$mz, .data$rt, .data$rtmin, .data$rtmax, .data$into,
            .data$maxo, .data$sample_index, .data$name, .data$formula,
            .data$adduct, .data$class, .data$grade
        )

    new(
        "LipidSearchSource",
        name = "lipid-search",
        version = res$version,
        metadata = sample_metadata,
        peaks = peaks,
        compounds = compounds
    )
}
