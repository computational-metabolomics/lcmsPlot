#' Define the options to use when plotting LC-MS data coming from
#' Compound Discoverer results.
#'
#' @param compounds_query A `character` value indicating the expression
#' used to filter compounds from the Compound Discoverer results.
#' The expression is evaluated on the compound table and can reference
#' the following columns:
#' \describe{
#'   \item{name}{Compound name.}
#'   \item{formula}{Chemical formula of the compound.}
#'   \item{adduct}{Ion adduct (e.g. `[M+H]+`, `[M-H]-`).}
#'   \item{rt}{Retention time of the compound (in seconds).}
#'   \item{rtmin}{Minimum retention time of the compound peak.}
#'   \item{rtmax}{Maximum retention time of the compound peak.}
#'   \item{mz}{Mass-to-charge ratio (m/z) of the detected ion.}
#'   \item{maxo}{Maximum observed peak intensity.}
#'   \item{into}{Integrated peak area reported by Compound Discoverer.}
#'   \item{checked}{Whether the compound is checked in the Compounds table.
#'     Available only once compounds have actually been checked in the Compound
#'     Discoverer application: `.cdResult` files store no such column until
#'     then, and the scripting node exports none. Referencing `checked` in that
#'     case raises an error rather than matching nothing.}
#' }
#' @param rt_extend A `numeric` value indicating how much (in seconds)
#' the retention time window should be extended on each side of the
#' compound peak when extracting and plotting chromatograms.
#' @return A function that takes an `lcmsPlot` object and returns a modified
#' version with the specified Compound Discoverer options stored in
#' `options$compound_discoverer`. It is intended for use with the `+` operator,
#' which incrementally layers new data or visual components onto
#' the `lcmsPlot` object.
#' @export
#' @examples
#' ## A minimal stand-in for a Compound Discoverer Scripting Node export.
#' ## See [CompoundDiscovererNodeSource()] for the table layout.
#' node_dir <- tempfile("cd_node")
#' dir.create(node_dir)
#'
#' write_tab <- function(df, name) {
#'     path <- file.path(node_dir, name)
#'     utils::write.table(
#'         df, path, sep = "\t", row.names = FALSE, quote = FALSE)
#'     path
#' }
#'
#' node_args <- list(Tables = list(
#'     list(DataFile = write_tab(
#'         data.frame(
#'             "Compounds ID" = c(1, 2),
#'             "Name" = c("L-Proline", "L-Kynurenine"),
#'             "Formula" = c("C5 H9 N O2", "C10 H12 N2 O3"),
#'             check.names = FALSE),
#'         "compounds.txt")),
#'     list(DataFile = write_tab(
#'         data.frame(
#'             "Compounds per File ID" = c(1, 2),
#'             "StudyFileID" = c(1, 2),
#'             "Area" = c(24050, 18700),
#'             "Intensity" = c(24050, 18700),
#'             check.names = FALSE),
#'         "cpf.txt")),
#'     list(DataFile = write_tab(
#'         data.frame(
#'             "Features ID" = c(1, 2),
#'             "mz" = c(116.0704, 209.0918),
#'             "Ion" = c("[M+H]+1", "[M+H]+1"),
#'             "Area" = c(24050, 18700),
#'             "RT [min]" = c(7.10, 6.10),
#'             "Left RT [min]" = c(6.95, 5.95),
#'             "Right RT [min]" = c(7.25, 6.25),
#'             check.names = FALSE),
#'         "features.txt")),
#'     list(DataFile = write_tab(
#'         data.frame(
#'             "Compounds ID" = c(1, 2),
#'             "Compounds per File ID" = c(1, 2),
#'             check.names = FALSE),
#'         "link_cmp_cpf.txt")),
#'     list(DataFile = write_tab(
#'         data.frame(
#'             "Compounds per File ID" = c(1, 2),
#'             "Features ID" = c(1, 2),
#'             check.names = FALSE),
#'         "link_cpf_feat.txt")),
#'     list(DataFile = write_tab(
#'         data.frame(
#'             "StudyFileID" = c(1, 2),
#'             "File Name" = c("l-proline-MS1.mzML", "l-kynurenine-MS1.mzML"),
#'             check.names = FALSE),
#'         "study_files.txt"))))
#'
#' json_path <- file.path(node_dir, "node_args.json")
#' writeLines(jsonlite::toJSON(node_args, auto_unbox = TRUE), json_path)
#'
#' mzml_dir <- file.path(tempdir(), "cd-node-example")
#' utils::unzip(
#'     system.file("extdata", "standards-mzml.zip", package = "lcmsPlot"),
#'     exdir = mzml_dir)
#'
#' ds <- CompoundDiscovererNodeSource(
#'     node_args = json_path,
#'     sample_paths = file.path(
#'         mzml_dir, "mzml",
#'         c("l-proline-MS1.mzML", "l-kynurenine-MS1.mzML")))
#'
#' ## `compounds_query` is evaluated on the compound table.
#' p <- lcmsPlot(ds) +
#'   lp_compound_discoverer(
#'     compounds_query = 'name %in% c("L-Proline", "L-Kynurenine")',
#'     rt_extend = 15
#'   ) +
#'   lp_chromatogram(highlight_peaks = TRUE) +
#'   lp_grid(rows = "sample_id", cols = "name", free_x = TRUE) +
#'   lp_labels(title = "Compound Discoverer example", legend = "Sample") +
#'   lp_legend(position = "bottom")
#' p
lp_compound_discoverer <- function(compounds_query = NULL, rt_extend = 10) {
    args <- as.list(environment())

    make_interface_function(
        name = "lp_compound_discoverer",
        args_list = args,
        fn = function(obj) {
            .require_source(
                obj,
                function(source) {
                    is_cd_result(source) || is_cd_node_source(source)
                },
                "lp_compound_discoverer: The data object is not a Compound Discoverer results connection or scripting-node source")

            obj@options$compound_discoverer <- args

            return(obj)
        }
    )
}

#' Define the options to use when plotting LC-MS data coming from
#' a LipidSearch result file.
#'
#' @param lipids_query A `character` value indicating the expression used to
#' select which lipids to plot. The expression is evaluated on the lipid table
#' and can reference the following columns:
#' \describe{
#'   \item{name}{Plot label: the lipid ion, disambiguated with its retention
#'   time when the same ion is reported more than once.}
#'   \item{lipid_ion}{Lipid ion as reported by LipidSearch, e.g. `AEA(18:2)+H`
#'   (4.2 `LipidIon` / 5.2 `LipidID`).}
#'   \item{lipid_group}{LipidSearch lipid group.}
#'   \item{class}{Lipid class, e.g. `PC`, `AcCa`.}
#'   \item{sub_class}{Lipid subclass (5.2 only; `NA` for 4.2).}
#'   \item{fatty_acid}{Fatty acid composition (4.2 only; `NA` for 5.2).}
#'   \item{formula}{Ion formula.}
#'   \item{adduct}{Adduct: the explicit `AdductIon` column for 5.2 (e.g. `M+H`),
#'   otherwise derived from the ion name for 4.2 (e.g. `+H`).}
#'   \item{calc_mz}{Theoretical m/z.}
#'   \item{rej}{LipidSearch rejection flag as a `logical` (`TRUE` if rejected).}
#'   \item{mz}{m/z used for extraction (observed where available).}
#'   \item{rt, rtmin, rtmax}{Retention time and peak bounds, in seconds.}
#'   \item{into}{Integrated peak area; `NA` where nothing was detected.}
#'   \item{maxo}{Peak height; `NA` where nothing was detected.}
#'   \item{grade}{LipidSearch identification grade (`A` to `D`).}
#'   \item{mscore, sn}{LipidSearch m-score / ID score and signal-to-noise
#'   (`sn` is 4.2 only).}
#'   \item{detected}{Whether LipidSearch integrated a peak in this sample.}
#'   \item{lipid_rank}{Dense rank of lipids by descending total area, for
#'   top-N selection.}
#' }
#'
#' The query selects *lipids*, not lipid/sample rows: a predicate such as
#' `grade == "A"` keeps the lipids graded A in at least one sample, and each of
#' them is then extracted from **every** sample - including samples that the
#' result file does not cover, which use the lipid's consensus m/z and
#' retention-time window.
#' @param rt_extend A `numeric` value indicating how much (in seconds) the
#' retention time window should be extended on each side of the lipid peak when
#' extracting and plotting chromatograms.
#' @return A function that takes an `lcmsPlot` object and returns a modified
#' version with the specified LipidSearch options stored in
#' `options$lipid_search`. It is intended for use with the `+` operator,
#' which incrementally layers new data or visual components onto
#' the `lcmsPlot` object.
#' @seealso [LipidSearchSource()]
#' @export
#' @examples
#' ## A minimal LipidSearch 4.2 export; see [LipidSearchSource()] for the
#' ## format details.
#' results_path <- file.path(tempdir(), "lipids_pos.txt")
#'
#' columns <- c(
#'     "Rej.", "LipidIon", "LipidGroup", "Class", "FattyAcid", "CalcMz",
#'     "IonFormula",
#'     "Area[c-1]", "Area[c-2]", "Height[c-1]", "Height[c-2]",
#'     "Rt[c-1]", "Rt[c-2]", "ObsMz[c-1]", "ObsMz[c-2]",
#'     "Grade[c-1]", "Grade[c-2]")
#'
#' rows <- list(
#'     c(0, "PRO(0:0)+H", "PRO(0:0)", "PRO", "(0:0)", 116.0706,
#'       "C5 H10 O2 N1", 24050, 0, 24050, 0, 7.10, 0,
#'       "116.0704", "", "A", ""),
#'     c(0, "KYN(0:0)+H", "KYN(0:0)", "KYN", "(0:0)", 209.0921,
#'       "C10 H13 O3 N2", 0, 18700, 0, 18700, 0, 6.10,
#'       "", "209.0918", "", "A"))
#'
#' writeLines(
#'     c("#[c-1]:l-proline-MS1.raw",
#'       "#[c-2]:l-kynurenine-MS1.raw",
#'       "#mScoreThreshold:5.0",
#'       "",
#'       vapply(
#'           c(list(columns), rows),
#'           function(x) paste0(paste(x, collapse = "\t"), "\t"),
#'           character(1))),
#'     results_path)
#'
#' mzml_dir <- file.path(tempdir(), "lipid-search-example")
#' utils::unzip(
#'     system.file("extdata", "standards-mzml.zip", package = "lcmsPlot"),
#'     exdir = mzml_dir)
#'
#' ds <- LipidSearchSource(
#'     results_path = results_path,
#'     sample_paths = file.path(
#'         mzml_dir, "mzml",
#'         c("l-proline-MS1.mzML", "l-kynurenine-MS1.mzML")))
#'
#' ## `lipids_query` selects lipids, which are then extracted from every sample.
#' p <- lcmsPlot(ds) +
#'   lp_lipid_search(
#'     lipids_query = 'grade == "A"',
#'     rt_extend = 30
#'   ) +
#'   lp_chromatogram(highlight_peaks = TRUE) +
#'   lp_grid(rows = "sample_id", cols = "name", free_x = TRUE) +
#'   lp_labels(title = "LipidSearch example", legend = "Sample") +
#'   lp_legend(position = "bottom")
#' p
lp_lipid_search <- function(lipids_query = NULL, rt_extend = 30) {
    args <- as.list(environment())

    make_interface_function(
        name = "lp_lipid_search",
        args_list = args,
        fn = function(obj) {
            .require_source(
                obj,
                is_lipid_search_source,
                "lp_lipid_search: The data object is not a LipidSearch source")

            obj@options$lipid_search <- args

            return(obj)
        }
    )
}

#' Define the options to use when plotting LC-MS data coming from
#' an mzTab-M file.
#'
#' @param features_query A `character` value indicating the expression used to
#' select which features to plot. The expression is evaluated on the feature
#' table and can reference the following columns:
#' \describe{
#'   \item{feature_id}{The feature's `SMF_ID`.}
#'   \item{name}{Compound name (`chemical_name` of the SML row referencing the
#'   feature); `NA` for unannotated features.}
#'   \item{formula}{Chemical formula (`chemical_formula`).}
#'   \item{database_identifier}{Database identifier of the compound.}
#'   \item{reliability}{Identification reliability, e.g. `1` to `4`.}
#'   \item{sml_id}{The `SML_ID` of the annotation.}
#'   \item{ambiguous}{`TRUE` when several SML rows reference the feature.}
#'   \item{adduct_ion, charge}{Adduct ion and charge of the feature.}
#'   \item{mz}{Mass-to-charge ratio (`exp_mass_to_charge`).}
#'   \item{rt, rtmin, rtmax}{Retention time and its bounds, in seconds.}
#'   \item{into}{Abundance in the sample; `NA` where the feature has none.}
#'   \item{detected}{Whether the feature has an abundance in the sample.}
#' }
#'
#' The query selects *features*, not feature/sample rows: a predicate such as
#' `into > 1e6` keeps the features above that abundance in at least one
#' sample, and each of them is then extracted from **every** sample.
#' @param rt_extend A `numeric` value indicating how much (in seconds) the
#' retention time window should be extended on each side of the feature when
#' extracting and plotting chromatograms.
#' @return A function that takes an `lcmsPlot` object and returns a modified
#' version with the specified mzTab-M options stored in `options$mztab`. It is
#' intended for use with the `+` operator, which incrementally layers new data
#' or visual components onto the `lcmsPlot` object.
#' @seealso [MzTabMSource()]
#' @export
#' @examples
#' ## A minimal mzTab-M file; see [MzTabMSource()] for the format details.
#' mzml_dir <- file.path(tempdir(), "mztab-layer-example")
#' utils::unzip(
#'     system.file("extdata", "standards-mzml.zip", package = "lcmsPlot"),
#'     exdir = mzml_dir)
#'
#' mztab_path <- file.path(tempdir(), "standards-layer.mzTab")
#' writeLines(
#'     c(
#'         "MTD\tms_run[1]-location\tfile:///data/l-proline-MS1.mzML",
#'         "MTD\tms_run[2]-location\tfile:///data/l-kynurenine-MS1.mzML",
#'         "MTD\tassay[1]-ms_run_ref\tms_run[1]",
#'         "MTD\tassay[2]-ms_run_ref\tms_run[2]",
#'         paste(
#'             "SMH", "SML_ID", "SMF_ID_REFS", "chemical_name",
#'             "reliability", sep = "\t"),
#'         paste("SML", 1, 1, "L-Proline", 1, sep = "\t"),
#'         paste("SML", 2, 2, "L-Kynurenine", 3, sep = "\t"),
#'         paste(
#'             "SFH", "SMF_ID", "exp_mass_to_charge",
#'             "retention_time_in_seconds", "retention_time_in_seconds_start",
#'             "retention_time_in_seconds_end", "abundance_assay[1]",
#'             "abundance_assay[2]", sep = "\t"),
#'         paste("SMF", 1, 116.0706, 426, 417, 435, 24050, "null", sep = "\t"),
#'         paste("SMF", 2, 209.0921, 366, 357, 375, "null", 18700, sep = "\t")),
#'     mztab_path)
#'
#' ds <- MzTabMSource(
#'     mztab_path,
#'     sample_paths = file.path(
#'         mzml_dir, "mzml",
#'         c("l-proline-MS1.mzML", "l-kynurenine-MS1.mzML")))
#'
#' ## Only the confidently identified features.
#' p <- lcmsPlot(ds) +
#'   lp_mztab(features_query = "reliability <= 2", rt_extend = 30) +
#'   lp_chromatogram(highlight_peaks = TRUE) +
#'   lp_grid(rows = "sample_id", cols = "name", free_x = TRUE) +
#'   lp_labels(title = "mzTab-M example", legend = "Sample") +
#'   lp_legend(position = "bottom")
#' p
lp_mztab <- function(features_query = NULL, rt_extend = 10) {
    args <- as.list(environment())

    make_interface_function(
        name = "lp_mztab",
        args_list = args,
        fn = function(obj) {
            .require_source(
                obj,
                is_mztab_source,
                "lp_mztab: The data object is not an mzTab-M source")

            obj@options$mztab <- args

            return(obj)
        }
    )
}
