# Default column-name candidates for the Compound Discoverer scripting-node
# table exports. CD column names vary slightly between versions and requested
# column sets, so each logical field is resolved against a list of candidates.
.cd_node_cols <- list(
    compound_id = c("Compounds ID"),
    name = c("Name"),
    formula = c("Formula", "ElementalCompositionFormula",
                "Elemental Composition Formula"),
    # Compound-level m/z (consolidated compound), used by the Compound-centric
    # "Unknown Compounds" export where features carry no m/z of their own. This
    # mirrors ConsolidatedUnknownCompoundItems.MassOverCharge in the DB reader.
    compound_mz = c("m/z", "mz", "MassOverCharge"),
    cpf_id = c("Compounds per File ID"),
    study_file_id = c("StudyFileID", "Study File ID"),
    # "Area Ref Ion" is the reference-ion area (the ion actually plotted);
    # "Area All Ions" is the summed area. "Area Max" appears on compound/feature
    # tables. Plain "Area" is kept first for older per-file-feature exports.
    cpf_area = c("Area", "Area Ref Ion", "Area All Ions", "Area Max"),
    cpf_intensity = c("Intensity", "Height", "Max. Intensity", "Intensity Max"),
    feature_id = c("Features ID"),
    feature_mz = c("mz", "m/z", "MassOverCharge"),
    feature_ion = c("Ion"),
    feature_area = c("Area", "Area Max"),
    # Adduct label as reported directly on the compound / instance tables.
    compound_adduct = c("Reference Ion"),
    rt = c("RT [min]", "Apex RT [min]", "RT in min", "RT", "Retention Time",
           "Ion_RT"),
    rt_left = c("Left RT [min]", "LeftRT", "RT left [min]", "Left RT"),
    rt_right = c("Right RT [min]", "RightRT", "RT right [min]", "Right RT"),
    # Peak full-width-at-half-maximum, used to derive the RT window when explicit
    # Left/Right RT bounds are not exported (CD 3.5 "Unknown Compounds").
    fwhm = c("FWHM in min", "FWHM [min]", "FWHM"),
    file_name = c("File Name", "FileName", "Physical File Name",
                  "Spectrum File Name", "Study File")
)

# Adduct strings (in priority order) used to select a representative
# "molecular" ion for a compound, matching the logic in CompoundsMZ.R.
.cd_node_molecular_ions <- c("[M+H]+1", "[M-H]-1")

#' Read the Compound Discoverer scripting-node tables
#'
#' Parses the `node_args.json` file produced by a Compound Discoverer
#' Scripting Node (or accepts an already-parsed list) and reads each referenced
#' tab-delimited table export into a `tibble`.
#'
#' @param node_args A `character` path to a `node_args.json` file or an
#' already-parsed `list` with the same structure.
#' @return A `list` of `tibble`s, one per entry in `Tables`.
#' @keywords internal
.cd_node_read_tables <- function(node_args) {
    if (is.character(node_args)) {
        if (length(node_args) != 1 || !file.exists(node_args)) {
            stop("node_args must be a path to an existing node_args.json file.")
        }
        node_args <- jsonlite::fromJSON(node_args, simplifyVector = FALSE)
    }

    if (is.null(node_args$Tables) || length(node_args$Tables) == 0) {
        stop("node_args does not contain any 'Tables'.")
    }

    lapply(node_args$Tables, function(tbl_spec) {
        data_file <- tbl_spec$DataFile
        if (is.null(data_file) || !file.exists(data_file)) {
            stop("Table data file not found: ", data_file)
        }
        utils::read.table(
            data_file,
            header = TRUE,
            sep = "\t",
            check.names = FALSE,
            stringsAsFactors = FALSE
        ) |> as_tibble()
    })
}

#' Classify the Compound Discoverer scripting-node tables
#'
#' Identifies the Compounds, Compounds per File, Features, the two link tables,
#' and (optionally) a study-file table from a list of exported tables, based on
#' their columns rather than their position.
#'
#' @param tables A `list` of `tibble`s as returned by `.cd_node_read_tables()`.
#' @return A named `list` with elements `compounds`, `cpf`, `features`,
#' `link_cmp_cpf`, `link_cpf_feat`, `link_cmp_feat`, and `study_files` (any of
#' which may be `NULL`).
#'
#' @details Two export topologies are supported:
#' * **Per-file feature** (targeted-style): each `features` row carries `m/z`,
#'   `Ion` and RT, joined to `cpf` via a Compounds-per-File \eqn{\leftrightarrow}
#'   Features link (`link_cpf_feat`).
#' * **Compound-centric** (CD "Unknown Compounds"): the consolidated `compounds`
#'   table carries `m/z` (`compound_mz`), features/ions link to the *compound*
#'   (`link_cmp_feat`), and per-file RT/area live on the `cpf` instance table.
#'   In this topology `features`/`link_cpf_feat` are optional.
#' @keywords internal
.cd_node_classify_tables <- function(tables) {
    has_id <- function(df, field) {
        !is.na(first_matching_column(df, .cd_node_cols[[field]]))
    }

    classified <- list(
        compounds = NULL, cpf = NULL, features = NULL,
        link_cmp_cpf = NULL, link_cpf_feat = NULL, link_cmp_feat = NULL,
        study_files = NULL
    )

    for (df in tables) {
        is_cmp_id <- has_id(df, "compound_id")
        is_cpf_id <- has_id(df, "cpf_id")
        is_feat_id <- has_id(df, "feature_id")

        if (is.null(classified$link_cmp_cpf) && is_cmp_id && is_cpf_id &&
            !is_feat_id) {
            classified$link_cmp_cpf <- df
        } else if (is.null(classified$link_cpf_feat) && is_cpf_id &&
                   is_feat_id && !is_cmp_id) {
            classified$link_cpf_feat <- df
        } else if (is.null(classified$link_cmp_feat) && is_cmp_id &&
                   is_feat_id && !is_cpf_id) {
            # Compound <-> Features link (Compound-centric export).
            classified$link_cmp_feat <- df
        } else if (is.null(classified$features) && is_feat_id &&
                   has_id(df, "feature_ion")) {
            # A Features/ion table: identified by its ID + an Ion column. m/z is
            # optional here - in the Compound-centric export it lives on the
            # compound, and the ion table is used only for the adduct label.
            classified$features <- df
        } else if (is.null(classified$cpf) && is_cpf_id &&
                   has_id(df, "study_file_id")) {
            classified$cpf <- df
        } else if (is.null(classified$compounds) && is_cmp_id &&
                   has_id(df, "name")) {
            classified$compounds <- df
        } else if (is.null(classified$study_files) &&
                   has_id(df, "study_file_id") && has_id(df, "file_name") &&
                   !is_cmp_id && !is_feat_id) {
            classified$study_files <- df
        }
    }

    # `compounds`, `cpf` and their link are always required. `features` and the
    # Compounds-per-File <-> Features link are only required when the compound
    # table does not itself carry an m/z (i.e. the per-file-feature topology);
    # in the Compound-centric topology the m/z comes from the compound instead.
    has_compound_mz <- !is.null(classified$compounds) &&
        has_id(classified$compounds, "compound_mz")
    required <- c("compounds", "cpf", "link_cmp_cpf")
    if (!has_compound_mz) {
        required <- c(required, "features", "link_cpf_feat")
    }

    missing <- required[vapply(
        required, function(n) is.null(classified[[n]]), logical(1)
    )]
    if (length(missing) > 0) {
        seen <- vapply(seq_along(tables), function(i) {
            paste0("  [", i, "] ", paste(colnames(tables[[i]]), collapse = " | "))
        }, character(1))
        stop(
            "Could not identify required Compound Discoverer tables: ",
            paste(missing, collapse = ", "),
            ".\nTables seen (columns):\n", paste(seen, collapse = "\n"),
            "\nEnsure the Scripting Node requests the Compounds, Compounds per ",
            "File, and their link table; and either a Features table with m/z ",
            "and Ion (+ its link), or an m/z column on the Compounds table."
        )
    }

    classified
}

#' Build the consolidated per-compound, per-file compound table
#'
#' Associate each compound, in each study file, with a representative ion
#' (preferring the molecular ion, otherwise the most abundant feature).
#'
#' @param tabs A `list` of classified tables from `.cd_node_classify_tables()`.
#' @return A `tibble` with one row per compound and study file, containing
#' `compound_id` (the Compound Discoverer "Compounds ID" primary key),
#' `study_file_id`, `name`, `formula`, `adduct`, `mz`, `rt`, `rtmin`, `rtmax`,
#' `into`, `maxo`, and `compound_rank` (a dense rank of compounds by descending
#' total area, for top-N selection).
#' @details Dispatches on the export topology: when a per-file `features` table
#' and its Compounds-per-File link are present the feature-level builder is used;
#' otherwise the Compound-centric builder is used (m/z from the compound, per-file
#' RT/area from the instance table).
#' @keywords internal
.cd_node_build_compounds <- function(tabs) {
    rows <- if (!is.null(tabs$features) && !is.null(tabs$link_cpf_feat)) {
        .cd_node_build_compounds_feature(tabs)
    } else {
        .cd_node_build_compounds_compound(tabs)
    }
    add_compound_rank(rows)
}

#' Synthesize unique, non-empty compound labels
#'
#' Compound Discoverer "Unknown Compounds" have an empty `Name`; the plot keys
#' compound identity on `name`, so a stable, unique, non-empty label is required.
#' Uses `Name` when present, else `Formula`, else `Compound <id>`, and appends
#' `" [#<id>]"` only where a label would otherwise collide across compound IDs.
#'
#' @keywords internal
.cd_node_compound_labels <- function(ids, names, formulas) {
    nz <- function(x) !is.na(x) & nzchar(trimws(x))
    base <- ifelse(nz(names), names,
        ifelse(nz(formulas), formulas, paste0("Compound ", ids)))
    dup <- base %in% base[duplicated(base)]
    ifelse(dup, paste0(base, " [#", ids, "]"), base)
}

#' Feature-level compound builder (per-file feature topology)
#' @keywords internal
.cd_node_build_compounds_feature <- function(tabs) {
    col <- function(df, field) first_matching_column(df, .cd_node_cols[[field]])

    cmp <- tabs$compounds
    cpf <- tabs$cpf
    feats <- tabs$features
    l_cc <- tabs$link_cmp_cpf
    l_cf <- tabs$link_cpf_feat

    cmp_id_c <- col(cmp, "compound_id")
    name_c <- col(cmp, "name")
    formula_c <- col(cmp, "formula")

    cpf_id_c <- col(cpf, "cpf_id")
    sfid_c <- col(cpf, "study_file_id")
    cpf_area_c <- col(cpf, "cpf_area")
    cpf_int_c <- col(cpf, "cpf_intensity")
    cpf_rt_c <- col(cpf, "rt")

    feat_id_c <- col(feats, "feature_id")
    feat_mz_c <- col(feats, "feature_mz")
    feat_ion_c <- col(feats, "feature_ion")
    feat_area_c <- col(feats, "feature_area")
    feat_rt_c <- col(feats, "rt")
    feat_lrt_c <- col(feats, "rt_left")
    feat_rrt_c <- col(feats, "rt_right")

    rows <- list()

    for (ci in seq_len(nrow(cmp))) {
        cid <- cmp[[cmp_id_c]][ci]
        cmp_name <- cmp[[name_c]][ci]
        cmp_formula <- if (!is.na(formula_c)) cmp[[formula_c]][ci] else NA_character_

        cpf_ids <- l_cc[[cpf_id_c]][l_cc[[col(l_cc, "compound_id")]] == cid]
        cpf_rows <- cpf[cpf[[cpf_id_c]] %in% cpf_ids, , drop = FALSE]

        for (ri in seq_len(nrow(cpf_rows))) {
            this_cpf_id <- cpf_rows[[cpf_id_c]][ri]
            sfid <- cpf_rows[[sfid_c]][ri]
            area <- if (!is.na(cpf_area_c)) cpf_rows[[cpf_area_c]][ri] else NA_real_
            intensity <- if (!is.na(cpf_int_c)) {
                cpf_rows[[cpf_int_c]][ri]
            } else {
                NA_real_
            }

            feat_ids <- l_cf[[feat_id_c]][
                l_cf[[col(l_cf, "cpf_id")]] == this_cpf_id]
            f <- feats[feats[[feat_id_c]] %in% feat_ids, , drop = FALSE]
            if (nrow(f) == 0) next

            # Prefer a molecular ion, otherwise the most abundant feature.
            sel <- integer(0)
            for (ion in .cd_node_molecular_ions) {
                hits <- which(f[[feat_ion_c]] == ion)
                if (length(hits) > 0) {
                    sel <- hits
                    break
                }
            }
            if (length(sel) == 0) {
                sel <- if (!is.na(feat_area_c)) which.max(f[[feat_area_c]]) else 1L
            } else if (length(sel) > 1 && !is.na(feat_area_c)) {
                sel <- sel[which.max(f[[feat_area_c]][sel])]
            } else {
                sel <- sel[1]
            }

            mz <- f[[feat_mz_c]][sel]
            adduct <- f[[feat_ion_c]][sel]

            rt <- if (!is.na(feat_rt_c)) {
                f[[feat_rt_c]][sel]
            } else if (!is.na(cpf_rt_c)) {
                cpf_rows[[cpf_rt_c]][ri]
            } else {
                NA_real_
            }
            rtmin <- if (!is.na(feat_lrt_c)) f[[feat_lrt_c]][sel] else rt
            rtmax <- if (!is.na(feat_rrt_c)) f[[feat_rrt_c]][sel] else rt

            rows[[length(rows) + 1]] <- tibble(
                compound_id = cid,
                study_file_id = sfid,
                name = cmp_name,
                formula = cmp_formula,
                adduct = adduct,
                mz = as.numeric(mz),
                rt = as.numeric(rt),
                rtmin = as.numeric(rtmin),
                rtmax = as.numeric(rtmax),
                into = as.numeric(area),
                maxo = as.numeric(intensity)
            )
        }
    }

    if (length(rows) == 0) {
        stop("No compound/feature associations were found in the export.")
    }

    do.call(rbind, rows)
}

#' Compound-centric compound builder (CD "Unknown Compounds" topology)
#'
#' Builds the per-compound, per-file table when m/z lives on the consolidated
#' `compounds` table and per-file RT/area live on the `cpf` (instance) table. The
#' representative adduct label is taken from the molecular ion of the compound's
#' `features`/ions via `link_cmp_feat` when available, otherwise `NA`.
#' @keywords internal
.cd_node_build_compounds_compound <- function(tabs) {
    col <- function(df, field) first_matching_column(df, .cd_node_cols[[field]])

    cmp <- tabs$compounds
    cpf <- tabs$cpf
    feats <- tabs$features        # optional
    l_cc <- tabs$link_cmp_cpf
    l_cmf <- tabs$link_cmp_feat   # optional

    cmp_id_c <- col(cmp, "compound_id")
    name_c <- col(cmp, "name")
    formula_c <- col(cmp, "formula")
    cmp_mz_c <- col(cmp, "compound_mz")
    cmp_adduct_c <- col(cmp, "compound_adduct")

    cpf_id_c <- col(cpf, "cpf_id")
    sfid_c <- col(cpf, "study_file_id")
    cpf_area_c <- col(cpf, "cpf_area")
    cpf_int_c <- col(cpf, "cpf_intensity")
    cpf_rt_c <- col(cpf, "rt")
    cpf_lrt_c <- col(cpf, "rt_left")
    cpf_rrt_c <- col(cpf, "rt_right")
    cpf_fwhm_c <- col(cpf, "fwhm")

    lcc_cmp_c <- col(l_cc, "compound_id")
    lcc_cpf_c <- col(l_cc, "cpf_id")

    # Optional adduct lookup: prefer the compound's own "Reference Ion", else the
    # (molecular) ion label from the linked features table.
    adduct_for <- function(cid) NA_character_
    if (!is.na(cmp_adduct_c)) {
        adduct_by_id <- setNames(cmp[[cmp_adduct_c]], as.character(cmp[[cmp_id_c]]))
        adduct_for <- function(cid) {
            a <- adduct_by_id[[as.character(cid)]]
            if (is.null(a) || is.na(a) || !nzchar(a)) NA_character_ else a
        }
    } else if (!is.null(feats) && !is.null(l_cmf)) {
        feat_id_c <- col(feats, "feature_id")
        feat_ion_c <- col(feats, "feature_ion")
        feat_area_c <- col(feats, "feature_area")
        lcmf_cmp_c <- col(l_cmf, "compound_id")
        lcmf_feat_c <- col(l_cmf, "feature_id")
        if (!is.na(feat_ion_c)) {
            adduct_for <- function(cid) {
                fids <- l_cmf[[lcmf_feat_c]][l_cmf[[lcmf_cmp_c]] == cid]
                f <- feats[feats[[feat_id_c]] %in% fids, , drop = FALSE]
                if (nrow(f) == 0) return(NA_character_)
                for (ion in .cd_node_molecular_ions) {
                    if (any(f[[feat_ion_c]] == ion, na.rm = TRUE)) return(ion)
                }
                if (!is.na(feat_area_c)) {
                    f[[feat_ion_c]][which.max(f[[feat_area_c]])]
                } else {
                    f[[feat_ion_c]][1]
                }
            }
        }
    }

    labels <- .cd_node_compound_labels(
        ids = cmp[[cmp_id_c]],
        names = if (!is.na(name_c)) cmp[[name_c]] else {
            rep(NA_character_, nrow(cmp))
        },
        formulas = if (!is.na(formula_c)) cmp[[formula_c]] else {
            rep(NA_character_, nrow(cmp))
        }
    )

    rows <- list()

    for (ci in seq_len(nrow(cmp))) {
        cid <- cmp[[cmp_id_c]][ci]
        label <- labels[ci]
        cmp_formula <- if (!is.na(formula_c)) cmp[[formula_c]][ci] else NA_character_
        mz <- if (!is.na(cmp_mz_c)) cmp[[cmp_mz_c]][ci] else NA_real_
        adduct <- adduct_for(cid)

        cpf_ids <- l_cc[[lcc_cpf_c]][l_cc[[lcc_cmp_c]] == cid]
        cpf_rows <- cpf[cpf[[cpf_id_c]] %in% cpf_ids, , drop = FALSE]
        if (nrow(cpf_rows) == 0) next

        for (ri in seq_len(nrow(cpf_rows))) {
            sfid <- cpf_rows[[sfid_c]][ri]
            rt <- if (!is.na(cpf_rt_c)) cpf_rows[[cpf_rt_c]][ri] else NA_real_
            # Prefer explicit Left/Right RT bounds; otherwise derive the peak
            # window from the FWHM around the apex; otherwise fall back to a
            # single RT point (rt_extend later widens it for extraction).
            if (!is.na(cpf_lrt_c)) {
                rtmin <- cpf_rows[[cpf_lrt_c]][ri]
                rtmax <- if (!is.na(cpf_rrt_c)) cpf_rows[[cpf_rrt_c]][ri] else rt
            } else if (!is.na(cpf_fwhm_c)) {
                fwhm <- cpf_rows[[cpf_fwhm_c]][ri]
                rtmin <- rt - fwhm / 2
                rtmax <- rt + fwhm / 2
            } else {
                rtmin <- rt
                rtmax <- rt
            }
            area <- if (!is.na(cpf_area_c)) cpf_rows[[cpf_area_c]][ri] else NA_real_
            intensity <- if (!is.na(cpf_int_c)) {
                cpf_rows[[cpf_int_c]][ri]
            } else {
                NA_real_
            }

            rows[[length(rows) + 1]] <- tibble(
                compound_id = cid,
                study_file_id = sfid,
                name = label,
                formula = cmp_formula,
                adduct = adduct,
                mz = as.numeric(mz),
                rt = as.numeric(rt),
                rtmin = as.numeric(rtmin),
                rtmax = as.numeric(rtmax),
                into = as.numeric(area),
                maxo = as.numeric(intensity)
            )
        }
    }

    if (length(rows) == 0) {
        stop("No compound/instance associations were found in the export.")
    }

    do.call(rbind, rows)
}

#' Map study files to sample paths
#'
#' Builds a sample metadata table from the distinct study files referenced in
#' the consolidated compound table and the user-supplied `sample_paths`. When a
#' study-file table with file names is available, study files are matched to
#' `sample_paths` by basename (ignoring extension);
#' otherwise study files are matched positionally by sorted study-file ID.
#'
#' @param compounds A `tibble` from `.cd_node_build_compounds()`.
#' @param study_files A `tibble` study-file table or `NULL`.
#' @param sample_paths A `character` vector of raw sample (e.g. mzML) paths.
#' @return A `tibble` with `study_file_id`, `sample_index`, `sample_id`, and
#' `sample_path`.
#' @keywords internal
.cd_node_build_metadata <- function(compounds, study_files, sample_paths) {
    sfids <- sort(unique(compounds$study_file_id))
    sample_basenames <- tools::file_path_sans_ext(basename(sample_paths))

    if (!is.null(study_files)) {
        sfid_c <- first_matching_column(study_files, .cd_node_cols$study_file_id)
        fname_c <- first_matching_column(study_files, .cd_node_cols$file_name)

        matched_path <- vapply(sfids, function(sfid) {
            fname <- study_files[[fname_c]][study_files[[sfid_c]] == sfid]
            if (length(fname) == 0) {
                stop("Study file ID ", sfid, " not found in study-file table.")
            }
            key <- tools::file_path_sans_ext(basename(fname[[1]]))
            idx <- match(key, sample_basenames)
            if (is.na(idx)) {
                stop(
                    "Could not match study file '", fname[[1]],
                    "' to any of the provided sample_paths."
                )
            }
            sample_paths[[idx]]
        }, character(1))
    } else {
        if (length(sfids) != length(sample_paths)) {
            stop(
                "No study-file table was found to match samples by name, and ",
                "the number of study files (", length(sfids), ") does not ",
                "equal the number of sample_paths (", length(sample_paths),
                "). Provide one sample path per study file, ordered by ",
                "study-file ID."
            )
        }
        matched_path <- sample_paths
    }

    tibble(
        study_file_id = sfids,
        sample_index = seq_along(sfids),
        sample_id = tools::file_path_sans_ext(basename(matched_path)),
        sample_path = matched_path
    )
}

#' Compound Discoverer scripting-node data source
#'
#' An S4 class wrapping the compound and peak data exported by a Compound
#' Discoverer Scripting Node.
#'
#' @slot name A `character` value identifying the data source.
#' @slot metadata A `data.frame` of sample metadata (one row per study file),
#' including a `sample_path` column with the raw (e.g. mzML) file paths.
#' @slot peaks A `data.frame` of detected peaks (one row per compound and study
#' file) with the standard `mz`, `rt`, `rtmin`, `rtmax`, `into`, `maxo`, and
#' `sample_index` columns, plus `name`, `formula`, and `adduct` annotations.
#' @slot compounds A `data.frame` of the consolidated per-compound, per-file
#' table used to drive chromatogram extraction.
#' @export
setClass(
    "CompoundDiscovererNodeSource",
    slots = list(
        name = "character",
        metadata = "data.frame",
        peaks = "data.frame",
        compounds = "data.frame"
    )
)

#' Show a summary of a `CompoundDiscovererNodeSource` object
#'
#' @param object An instance of class `CompoundDiscovererNodeSource`.
#' @return Invisible \code{NULL}
#' @export
setMethod(
    f = "show",
    signature = "CompoundDiscovererNodeSource",
    function(object) {
        cat("Object of class", class(object), "\n")
        cat(" Name:", object@name, "\n")
        cat(
            " Samples:", nrow(object@metadata), "\n",
            " Compounds:", length(unique(object@compounds$name)), "\n",
            " Peaks:", nrow(object@peaks), "rows\n"
        )
    }
)

#' Create a Compound Discoverer data source from a Scripting Node export
#'
#' Reads the `node_args.json` file produced by a Compound Discoverer Scripting
#' Node together with the tab-delimited table exports it references, and builds
#' a data source that can be plotted with `lcmsPlot()`. Chromatograms are later
#' extracted from the raw (mzML) files supplied via `sample_paths`.
#'
#' The Scripting Node must be configured to request the Compounds, Compounds per
#' File, and Features tables (with their link tables). For each compound and
#' study file a representative ion is selected, preferring a molecular ion
#' and otherwise the most abundant feature.
#'
#' @param node_args A `character` path to the `node_args.json` file passed by
#' Compound Discoverer or an already-parsed `list`.
#' @param sample_paths A `character` vector of paths to the raw sample files
#' (e.g. mzML) used to extract chromatograms. Study files are matched to these
#' paths by file basename when a study-file table is available, otherwise
#' positionally by study-file ID.
#' @param metadata A `data.frame` (optional) of additional sample metadata. It
#' must contain one row per sample, aligned with `sample_paths`. When supplied,
#' its columns are column-bound onto the generated metadata.
#' @return An object of class `CompoundDiscovererNodeSource`.
#' @seealso [lp_compound_discoverer()]
#' @export
#' @examples
#' \dontrun{
#' ## Inside a Compound Discoverer Scripting Node, the JSON path is the
#' ## 6th command-line argument:
#' node_args <- commandArgs()[6]
#'
#' ds <- CompoundDiscovererNodeSource(
#'     node_args = node_args,
#'     sample_paths = c("path/to/sample1.mzML", "path/to/sample2.mzML")
#' )
#'
#' lcmsPlot(ds) +
#'     lp_compound_discoverer(
#'         compounds_query = 'name %in% c("Betaine", "Creatine")',
#'         rt_extend = 5
#'     ) +
#'     lp_chromatogram(highlight_peaks = TRUE) +
#'     lp_grid(rows = "sample_id", cols = "name", free_x = TRUE)
#' }
CompoundDiscovererNodeSource <- function(
    node_args,
    sample_paths,
    metadata = NULL
) {
    tables <- .cd_node_read_tables(node_args)
    tabs <- .cd_node_classify_tables(tables)

    # Compound Discoverer reports retention times in minutes
    compounds <- .cd_node_build_compounds(tabs) |>
        convert_rt_to_seconds()

    sample_metadata <- .cd_node_build_metadata(
        compounds, tabs$study_files, sample_paths)

    if (!is.null(metadata)) {
        if (nrow(metadata) != length(sample_paths)) {
            stop("metadata must have one row per sample path.")
        }
        path_idx <- match(sample_metadata$sample_path, sample_paths)
        sample_metadata <- as_tibble(
            cbind(sample_metadata, as_tibble(metadata)[path_idx, , drop = FALSE])
        )
    }

    # Attach sample_index to the consolidated compound table and peaks.
    compounds <- compounds |>
        left_join(
            sample_metadata |> select(.data$study_file_id, .data$sample_index),
            by = "study_file_id"
        )

    peaks <- compounds |>
        select(
            .data$mz, .data$rt, .data$rtmin, .data$rtmax,
            .data$into, .data$maxo, .data$sample_index,
            .data$name, .data$formula, .data$adduct
        )

    new(
        "CompoundDiscovererNodeSource",
        name = "compound-discoverer-node",
        metadata = sample_metadata,
        peaks = peaks,
        compounds = compounds
    )
}

#' Write a Compound Discoverer `node_response.json` adding a per-compound column
#'
#' Writes a `node_response.json` that returns **only** the Compounds table (the
#' one table modified) with one new column appended, matched on the Compounds
#' table's `"Compounds ID"` primary key. Each row's value comes from
#' `id_to_value` (rows without a mapping get an empty string). The column is
#' displayed in Compound Discoverer via the supplied `SpecialCellRenderer` GUID -
#' e.g. CD's filename renderer, which turns the cell into a clickable file link.
#' Only the modified table is returned.
#'
#' @param node_args A `character` path to the `node_args.json` file.
#' @param id_to_value A `data.frame` with columns `compound_id` and `value`
#' mapping each Compounds ID to its cell value (e.g. a plot file path).
#' @param column_name The new column's header. Default `"Plot"`.
#' @param renderer The `SpecialCellRenderer` GUID string.
#' @param position_after An existing Compounds column to place the new column
#' after. Default `"Name"`.
#' @return The path to the written `node_response.json` (invisibly), or `NULL` if
#' the `node_args` has no `ExpectedResponsePath` (e.g. a standalone run).
#' @keywords internal
.cd_node_write_plot_column <- function(
    node_args,
    id_to_value,
    column_name = "Plot",
    renderer = "EB29D794-4F2E-4785-8B80-A24D8C0FB3E4",
    position_after = "Name"
) {
    if (!is.character(node_args) || length(node_args) != 1 ||
        !file.exists(node_args)) {
        stop("node_args must be a path to an existing node_args.json file.")
    }
    spec <- jsonlite::fromJSON(node_args, simplifyVector = FALSE)

    response_path <- spec$ExpectedResponsePath
    if (is.null(response_path) || !nzchar(response_path)) {
        return(invisible(NULL))
    }

    # Locate the Compounds table entry.
    table_names <- vapply(spec$Tables, function(t) {
        if (is.null(t$TableName)) "" else t$TableName
    }, character(1))
    cmp_idx <- match("Compounds", table_names)
    if (is.na(cmp_idx)) {
        stop("node_args does not contain a 'Compounds' table.")
    }
    cmp_tbl <- spec$Tables[[cmp_idx]]
    data_file <- cmp_tbl$DataFile
    if (is.null(data_file) || !file.exists(data_file)) {
        stop("Compounds table data file not found: ", data_file)
    }

    # Read the export as text so the original values are preserved verbatim
    # (CD quotes every field, including numbers, so read all as character and
    # rewrite quoted to reproduce the exact format without reformatting numbers).
    df <- utils::read.table(
        data_file, header = TRUE, sep = "\t", check.names = FALSE,
        stringsAsFactors = FALSE, colClasses = "character",
        quote = "\"", comment.char = ""
    )

    id_col <- first_matching_column(df, .cd_node_cols$compound_id)
    if (is.na(id_col)) {
        stop("Compounds export has no recognised 'Compounds ID' column.")
    }

    lut <- stats::setNames(
        as.character(id_to_value$value),
        as.character(id_to_value$compound_id)
    )
    vals <- unname(lut[as.character(df[[id_col]])])
    vals[is.na(vals)] <- ""
    df[[column_name]] <- vals

    out_file <- sub("\\.txt$", ".out.txt", data_file)
    if (identical(out_file, data_file)) {
        out_file <- paste0(data_file, ".out.txt")
    }
    utils::write.table(
        df, out_file, sep = "\t", row.names = FALSE, quote = TRUE)

    # Append the column description and repoint the Compounds DataFile.
    new_col <- list(
        ColumnName = column_name,
        ID = "",
        DataType = "String",
        Options = list(
            PositionAfter = position_after,
            SpecialCellRenderer = renderer
        )
    )
    cmp_tbl$ColumnDescriptions[[length(cmp_tbl$ColumnDescriptions) + 1]] <-
        new_col
    cmp_tbl$DataFile <- out_file
    # Return only the modified Compounds table.
    spec$Tables <- list(cmp_tbl)

    json <- jsonlite::toJSON(
        spec, auto_unbox = TRUE, pretty = TRUE, null = "null")
    # jsonlite serializes empty objects ({}) as empty arrays ([]); restore the
    # object form for the fields Compound Discoverer expects to be objects.
    json <- gsub(
        '("Options"|"NodeParameters")(\\s*):(\\s*)\\[\\s*\\]',
        "\\1\\2:\\3{}", json, perl = TRUE
    )

    writeLines(json, response_path)
    invisible(response_path)
}
