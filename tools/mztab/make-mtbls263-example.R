# Build the mzTab-M example shipped in inst/extdata from MetaboLights MTBLS263.
#
# The HUPO-PSI mzTab-M example MTBLS263.mztab (Progenesis QI results) refers to
# the six positive-mode runs of MetaboLights study MTBLS263. This script keeps
# the features eluting between 380 and 680 seconds, and crops each run to their
# retention time bounds (widened by `rt_margin`) and m/z, so that the example
# stays small.
#
# Run from the package root:
#   Rscript tools/mztab/make-mtbls263-example.R <mztab> <raw_dir>
#
# where <mztab> is examples/2_0-Metabolomics_Release/MTBLS263.mztab from
# https://github.com/HUPO-PSI/mzTab, and <raw_dir> holds the *_POS.mzML files
# from https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS263/FILES/

args <- commandArgs(trailingOnly = TRUE)
mztab_path <- args[[1]]
raw_dir <- args[[2]]

rt_range <- c(380, 680)
rt_margin <- 60
mz_tolerance <- 0.05
out_dir <- file.path("inst", "extdata")

# Features, with the summary and evidence rows that refer to them.
lines <- readLines(mztab_path)
fields <- strsplit(lines, "\t", fixed = TRUE)
prefix <- vapply(fields, `[[`, character(1), 1)
field <- function(i, k) fields[[i]][[k]]

smf_header <- fields[[which(prefix == "SFH")]]
col <- function(name) match(name, smf_header)
smf_rows <- which(prefix == "SMF")
smf_rt <- as.numeric(vapply(
    smf_rows, field, character(1), col("retention_time_in_seconds")))
kept_smf <- smf_rows[smf_rt >= rt_range[1] & smf_rt <= rt_range[2]]
smf_ids <- vapply(kept_smf, field, character(1), col("SMF_ID"))
feature_mz <- as.numeric(vapply(
    kept_smf, field, character(1), col("exp_mass_to_charge")))
crop_rt <- c(
    min(as.numeric(vapply(
        kept_smf, field, character(1),
        col("retention_time_in_seconds_start")))) - rt_margin,
    max(as.numeric(vapply(
        kept_smf, field, character(1),
        col("retention_time_in_seconds_end")))) + rt_margin)
sme_ids <- unlist(strsplit(
    vapply(kept_smf, field, character(1), col("SME_ID_REFS")), "\\s*\\|\\s*"))

refs <- function(rows, k) strsplit(
    vapply(rows, field, character(1), k), "\\s*\\|\\s*")
sml_rows <- which(prefix == "SML")
kept_sml <- sml_rows[vapply(
    refs(sml_rows, 3), function(r) any(r %in% smf_ids), logical(1))]
sme_rows <- which(prefix == "SME")
kept_sme <- sme_rows[vapply(sme_rows, field, character(1), 2) %in% sme_ids]

drop <- setdiff(c(sml_rows, smf_rows, sme_rows), c(kept_sml, kept_smf, kept_sme))
writeLines(lines[-drop], file.path(out_dir, "MTBLS263.mztab"))

# The raw data of those features only.
mzml_dir <- file.path(tempfile("mtbls263"), "mzml")
dir.create(mzml_dir, recursive = TRUE)
for (path in list.files(raw_dir, "_POS\\.mzML$", full.names = TRUE)) {
    message("Cropping ", basename(path))
    sps <- Spectra::Spectra(path, backend = Spectra::MsBackendMzR())
    sps <- Spectra::filterMsLevel(sps, 1L)
    sps <- Spectra::filterRt(sps, crop_rt)
    sps <- Spectra::filterMzValues(
        sps, sort(feature_mz), tolerance = mz_tolerance)
    sps <- Spectra::setBackend(sps, Spectra::MsBackendMemory())
    Spectra::export(
        sps,
        Spectra::MsBackendMzR(),
        file = file.path(mzml_dir, basename(path)),
        copy = TRUE)
}

zip_path <- file.path(normalizePath(out_dir), "mtbls263-mzml.zip")
unlink(zip_path)
withr::with_dir(dirname(mzml_dir), utils::zip(
    zip_path, file.path("mzml", list.files(mzml_dir)), flags = "-9Xq"))
