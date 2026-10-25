#' Read the sections of an mzTab-M file
#'
#' mzTab-M 2.0 files are tab-separated, with the section of each line given by
#' its first field: `MTD` (metadata), `SMH`/`SML` (small molecule summary
#' header/rows), `SFH`/`SMF` (small molecule feature header/rows), `SEH`/`SME`
#' (evidence) and `COM` (comments). Only the metadata, summary and feature
#' sections are read.
#'
#' @param path A `character` path to an mzTab-M file.
#' @return A `list` with `mtd`, a named `character` vector of metadata values
#' keyed by metadata name (e.g. `ms_run[1]-location`), and `sml`/`smf`, the
#' summary and feature tables (`NULL` when the section is absent).
#' @keywords internal
.mztab_read <- function(path) {
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    prefix <- sub("\t.*$", "", lines)

    mtd_fields <- strsplit(lines[prefix == "MTD"], "\t", fixed = TRUE)
    mtd_fields <- mtd_fields[lengths(mtd_fields) >= 3]
    mtd <- vapply(mtd_fields, function(f) trimws(f[[3]]), character(1))
    names(mtd) <- vapply(mtd_fields, function(f) trimws(f[[2]]), character(1))

    list(
        mtd = mtd,
        sml = .mztab_section(lines, prefix, "SMH", "SML"),
        smf = .mztab_section(lines, prefix, "SFH", "SMF")
    )
}

#' Parse one table section of an mzTab-M file
#'
#' Rows are cut or padded to the header's width, as exports often end lines
#' with stray tabs. The literal `null` (mzTab's missing value) and empty fields
#' become `NA`, and each column is converted to its natural type.
#'
#' @param lines A `character` vector of the file's lines.
#' @param prefix A `character` vector of the first field of each line.
#' @param header_prefix The prefix of the section header line, e.g. `"SFH"`.
#' @param row_prefix The prefix of the section rows, e.g. `"SMF"`.
#' @return A `tibble`, or `NULL` if the file has no such section.
#' @keywords internal
.mztab_section <- function(lines, prefix, header_prefix, row_prefix) {
    header <- lines[prefix == header_prefix]
    if (length(header) == 0) {
        return(NULL)
    }

    columns <- trimws(strsplit(header[[1]], "\t", fixed = TRUE)[[1]])
    fields <- lapply(
        strsplit(lines[prefix == row_prefix], "\t", fixed = TRUE),
        function(f) {
            length(f) <- length(columns)
            f
        })
    values <- if (length(fields) > 0) {
        do.call(rbind, fields)
    } else {
        matrix(character(), nrow = 0, ncol = length(columns))
    }

    # The first column is the line prefix; unnamed columns come from trailing
    # tabs on the header line.
    keep <- nzchar(columns) & seq_along(columns) > 1
    df <- as.data.frame(values[, keep, drop = FALSE], stringsAsFactors = FALSE)
    colnames(df) <- columns[keep]

    df[] <- lapply(df, function(x) {
        x <- trimws(x)
        x[x %in% c("null", "")] <- NA
        utils::type.convert(x, as.is = TRUE)
    })

    as_tibble(df)
}

#' Map the assays of an mzTab-M file to the supplied sample paths
#'
#' Each `assay[n]-ms_run_ref` names an `ms_run[m]`, whose `-location` is the
#' URI of the raw file. The location's basename is matched to the basenames of
#' `sample_paths`, falling back to a match without the file extension (so a
#' result file listing `.raw` files works with converted `.mzML` files).
#'
#' @param mtd A named `character` vector of metadata values.
#' @param sample_paths A `character` vector of raw sample file paths.
#' @return A `tibble` with `assay_index`, `ms_run` and `sample_index` columns.
#' @keywords internal
.mztab_assay_runs <- function(mtd, sample_paths) {
    run_refs <- mtd[grepl("^assay\\[\\d+\\]-ms_run_ref$", names(mtd))]
    if (length(run_refs) == 0) {
        stop("The mzTab-M file declares no 'assay[n]-ms_run_ref' metadata.")
    }

    refs <- strsplit(run_refs, "\\s*\\|\\s*")
    multi_run <- lengths(refs) != 1
    if (any(multi_run)) {
        stop(
            "Assays referencing several ms_runs are not supported: ",
            paste(names(run_refs)[multi_run], collapse = ", "))
    }

    ms_run <- unlist(refs, use.names = FALSE)
    locations <- unname(mtd[paste0(ms_run, "-location")])
    if (anyNA(locations)) {
        stop(
            "The mzTab-M file declares no location for: ",
            paste(unique(ms_run[is.na(locations)]), collapse = ", "))
    }

    files <- vapply(
        locations,
        function(x) basename(gsub("\\\\", "/", utils::URLdecode(x))),
        character(1),
        USE.NAMES = FALSE)

    sample_files <- basename(sample_paths)
    sample_index <- match(files, sample_files)
    no_match <- is.na(sample_index)
    sample_index[no_match] <- match(
        tools::file_path_sans_ext(files[no_match]),
        tools::file_path_sans_ext(sample_files))

    if (anyNA(sample_index)) {
        stop(
            "Could not match these mzTab-M ms_run locations to sample_paths: ",
            paste(unique(files[is.na(sample_index)]), collapse = ", "))
    }

    tibble(
        assay_index = as.integer(
            sub("^assay\\[(\\d+)\\].*$", "\\1", names(run_refs))),
        ms_run = ms_run,
        sample_index = sample_index
    )
}

#' Build the per-sample feature abundances from the SMF section
#'
#' @param smf The SMF table.
#' @param assays The assay-to-sample map from `.mztab_assay_runs()`.
#' @return A `list` with `features` (one row per SMF feature, with its m/z and
#' retention time window) and `abundances` (one row per feature and sample).
#' @keywords internal
.mztab_features <- function(smf, assays) {
    if (is.null(smf)) {
        stop("The mzTab-M file has no small molecule feature (SFH/SMF) section.")
    }

    required <- c("SMF_ID", "exp_mass_to_charge", "retention_time_in_seconds")
    missing <- setdiff(required, colnames(smf))
    if (length(missing) > 0) {
        stop(
            "Missing required columns in the mzTab-M SMF section: ",
            paste(missing, collapse = ", "))
    }

    optional <- c(
        "adduct_ion", "charge",
        "retention_time_in_seconds_start", "retention_time_in_seconds_end")
    for (col in setdiff(optional, colnames(smf))) {
        smf[[col]] <- NA
    }

    abundance_cols <- grep(
        "^abundance_assay\\[\\d+\\]$", colnames(smf), value = TRUE)

    # Feature-level retention times, already in seconds. Without reported
    # bounds the window collapses to the apex and is widened by rt_extend.
    features <- smf |>
        transmute(
            feature_id = as.character(.data$SMF_ID),
            mz = as.numeric(.data$exp_mass_to_charge),
            rt = as.numeric(.data$retention_time_in_seconds),
            rtmin = dplyr::coalesce(
                as.numeric(.data$retention_time_in_seconds_start), .data$rt),
            rtmax = dplyr::coalesce(
                as.numeric(.data$retention_time_in_seconds_end), .data$rt),
            adduct_ion = as.character(.data$adduct_ion),
            charge = as.integer(.data$charge)
        )

    abundances <- smf |>
        transmute(
            feature_id = as.character(.data$SMF_ID),
            across(all_of(abundance_cols), as.numeric)
        ) |>
        tidyr::pivot_longer(
            all_of(abundance_cols),
            names_to = "assay_index",
            values_to = "into"
        ) |>
        mutate(assay_index = as.integer(
            sub("^abundance_assay\\[(\\d+)\\]$", "\\1", .data$assay_index))) |>
        filter(!is.na(.data$into)) |>
        inner_join(assays, by = "assay_index") |>
        select("feature_id", "sample_index", "into")

    list(features = features, abundances = abundances)
}

#' Build the feature annotations from the SML section
#'
#' Each SML row lists the features it summarises in `SMF_ID_REFS`
#' (`|`-separated). A feature referenced by several SML rows keeps the first,
#' and is flagged as `ambiguous`.
#'
#' @param sml The SML table, or `NULL`.
#' @return A `tibble` with one row per annotated feature: `feature_id`,
#' `sml_id`, `name`, `formula`, `database_identifier`, `reliability` and
#' `ambiguous`.
#' @keywords internal
.mztab_annotations <- function(sml) {
    if (is.null(sml) || !all(c("SML_ID", "SMF_ID_REFS") %in% colnames(sml))) {
        return(tibble(
            feature_id = character(),
            sml_id = character(),
            name = character(),
            formula = character(),
            database_identifier = character(),
            reliability = character(),
            ambiguous = logical()
        ))
    }

    optional <- c(
        "chemical_name", "chemical_formula", "database_identifier",
        "reliability")
    for (col in setdiff(optional, colnames(sml))) {
        sml[[col]] <- NA_character_
    }

    sml |>
        transmute(
            sml_id = as.character(.data$SML_ID),
            feature_id = strsplit(
                as.character(.data$SMF_ID_REFS), "\\s*\\|\\s*"),
            name = as.character(.data$chemical_name),
            formula = as.character(.data$chemical_formula),
            database_identifier = as.character(.data$database_identifier),
            reliability = .data$reliability
        ) |>
        tidyr::unnest("feature_id") |>
        filter(!is.na(.data$feature_id), nzchar(.data$feature_id)) |>
        group_by(.data$feature_id) |>
        mutate(ambiguous = n() > 1) |>
        ungroup() |>
        distinct(.data$feature_id, .keep_all = TRUE) |>
        relocate("feature_id")
}

#' mzTab-M data source
#'
#' An S4 class wrapping the small molecule features of an mzTab-M 2.0 file.
#' It extends `ExternalDataSource`, so its `peaks` hold the features detected
#' in each sample.
#'
#' @slot compounds A `data.frame` with one row per feature and sample -
#' including samples in which the feature has no abundance - used to drive
#' chromatogram extraction.
#' @seealso [MzTabMSource()], [lp_mztab()]
#' @export
setClass(
    "MzTabMSource",
    contains = "ExternalDataSource",
    slots = list(
        compounds = "data.frame"
    )
)

#' Show a summary of a `MzTabMSource` object
#'
#' @param object An instance of class `MzTabMSource`.
#' @return Invisible \code{NULL}
#' @export
setMethod(
    f = "show",
    signature = "MzTabMSource",
    function(object) {
        features <- object@compounds |>
            distinct(.data$feature_id, .keep_all = TRUE)

        cat("Object of class", class(object), "\n")
        cat(" Name:", object@name, "\n")
        cat(" Samples:", nrow(object@metadata), "\n")
        cat(
            " Features:", nrow(features),
            paste0("(", sum(!is.na(features$name)), " annotated)"), "\n")
        cat(" Peaks:", nrow(object@peaks), "rows\n")
    }
)

#' Create a data source from an mzTab-M file
#'
#' Reads the small molecule features (SMF) of an mzTab-M 2.0 file - as
#' exported by MS-DIAL, MZmine or OpenMS, among others - and pairs them with
#' the raw files they were computed from. Chromatograms are later extracted
#' from the raw files supplied via `sample_paths`, and the features can be
#' selected with [lp_mztab()].
#'
#' @section Sample-to-file mapping:
#' Each `assay[n]-ms_run_ref` of the metadata section names an `ms_run[m]`,
#' whose `ms_run[m]-location` is the URI of the raw file. The location's
#' basename is matched to the basenames of `sample_paths`, ignoring the file
#' extension when there is no exact match, so a file listing `.raw` files works
#' with converted `.mzML` files. Every run must match a sample path. Sample
#' paths that no run matches are still plotted, with no peaks to highlight.
#'
#' @section Feature columns:
#' Each SMF row is one feature across all assays. It gives one peak per assay
#' with an abundance:
#' \describe{
#'   \item{feature_id}{`SMF_ID`.}
#'   \item{mz}{`exp_mass_to_charge`.}
#'   \item{rt, rtmin, rtmax}{`retention_time_in_seconds` and its `_start`/`_end`
#'   bounds, shared by every sample. When the bounds are `null` they default to
#'   `rt`, so the extraction window is `rt ± rt_extend`.}
#'   \item{into}{`abundance_assay[n]`.}
#'   \item{maxo}{`NA`: mzTab-M reports no peak height.}
#'   \item{adduct_ion, charge}{As reported.}
#' }
#' Features are annotated from the small molecule summary (SML) rows that list
#' them in `SMF_ID_REFS`: `name` (`chemical_name`), `formula`
#' (`chemical_formula`), `database_identifier`, `reliability` and `sml_id`. A
#' feature referenced by several SML rows keeps the first, and its `ambiguous`
#' column is `TRUE`. Unannotated features have `NA` annotations.
#'
#' @param mztab_path A `character` path to an mzTab-M 2.0 file.
#' @param sample_paths A `character` vector of paths to raw sample files
#' (e.g., `mzML`).
#' @param metadata_path A `character` value (optional) indicating the
#' path to a metadata file (CSV or TSV), with one row per sample path.
#' If `NULL`, a metadata `data.frame` is created from `sample_paths`. A
#' `sample_id` column (the file name without extension) is added when the
#' metadata has none.
#' @return An object of class `MzTabMSource`, which extends
#' `ExternalDataSource`.
#' @seealso [lp_mztab()]
#' @export
#' @examples
#' mzml_dir <- file.path(tempdir(), "mztab-example")
#' utils::unzip(
#'     system.file("extdata", "standards-mzml.zip", package = "lcmsPlot"),
#'     exdir = mzml_dir)
#'
#' sample_paths <- file.path(
#'     mzml_dir, "mzml",
#'     c("l-proline-MS1.mzML", "l-kynurenine-MS1.mzML"))
#'
#' ## A minimal mzTab-M file: two runs, one feature found in each.
#' mztab_path <- file.path(tempdir(), "standards.mzTab")
#' writeLines(
#'     c(
#'         "MTD\tmzTab-version\t2.0.0-M",
#'         "MTD\tms_run[1]-location\tfile:///data/l-proline-MS1.mzML",
#'         "MTD\tms_run[2]-location\tfile:///data/l-kynurenine-MS1.mzML",
#'         "MTD\tassay[1]-ms_run_ref\tms_run[1]",
#'         "MTD\tassay[2]-ms_run_ref\tms_run[2]",
#'         "",
#'         paste(
#'             "SMH", "SML_ID", "SMF_ID_REFS", "chemical_name",
#'             "chemical_formula", "reliability", sep = "\t"),
#'         paste("SML", 1, 1, "L-Proline", "C5H9NO2", 1, sep = "\t"),
#'         paste("SML", 2, 2, "L-Kynurenine", "C10H12N2O3", 2, sep = "\t"),
#'         "",
#'         paste(
#'             "SFH", "SMF_ID", "adduct_ion", "exp_mass_to_charge", "charge",
#'             "retention_time_in_seconds", "retention_time_in_seconds_start",
#'             "retention_time_in_seconds_end", "abundance_assay[1]",
#'             "abundance_assay[2]", sep = "\t"),
#'         paste(
#'             "SMF", 1, "[M+H]1+", 116.0706, 1, 426, 417, 435, 24050, "null",
#'             sep = "\t"),
#'         paste(
#'             "SMF", 2, "[M+H]1+", 209.0921, 1, 366, 357, 375, "null", 18700,
#'             sep = "\t")),
#'     mztab_path)
#'
#' ds <- MzTabMSource(mztab_path, sample_paths)
#' ds
#'
#' p <- lcmsPlot(ds) +
#'     lp_mztab(rt_extend = 30) +
#'     lp_chromatogram(highlight_peaks = TRUE) +
#'     lp_grid(rows = "sample_id", cols = "name", free_x = TRUE)
#' p
MzTabMSource <- function(mztab_path, sample_paths, metadata_path = NULL) {
    metadata <- process_metadata(sample_paths, metadata_path)
    if (!"sample_id" %in% colnames(metadata)) {
        metadata$sample_id <- tools::file_path_sans_ext(basename(sample_paths))
    }

    res <- .mztab_read(mztab_path)
    assays <- .mztab_assay_runs(res$mtd, sample_paths)
    smf <- .mztab_features(res$smf, assays)
    annotations <- .mztab_annotations(res$sml)

    # The retention time window is feature-level, so every feature can be
    # extracted from every sample; `detected` marks the samples in which it
    # has an abundance.
    compounds <- smf$features |>
        cross_join(tibble(sample_index = seq_len(nrow(metadata)))) |>
        left_join(smf$abundances, by = c("feature_id", "sample_index")) |>
        left_join(annotations, by = "feature_id") |>
        mutate(
            ambiguous = dplyr::coalesce(.data$ambiguous, FALSE),
            maxo = NA_real_,
            detected = !is.na(.data$into)
        ) |>
        select(
            "feature_id", "name", "formula", "database_identifier",
            "reliability", "sml_id", "ambiguous", "adduct_ion", "charge",
            "mz", "rt", "rtmin", "rtmax", "into", "maxo", "detected",
            "sample_index"
        )

    peaks <- compounds |>
        filter(.data$detected) |>
        select(
            "mz", "rt", "rtmin", "rtmax", "into", "maxo", "sample_index",
            "feature_id", "name", "formula", "adduct_ion", "charge",
            "reliability"
        )

    new(
        "MzTabMSource",
        name = "mztab-m",
        metadata = metadata,
        peaks = peaks,
        compounds = compounds
    )
}
