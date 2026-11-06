# =============================================================================
# Install the R packages the Compounds Chromatograms node needs.
# -----------------------------------------------------------------------------
# Run by install.ps1 with the Rscript that Compound Discoverer will use:
#   Rscript install-packages.R <lcmsPlot version> [git ref]
#
# lcmsPlot is installed from the GitHub release tag v<version> that matches
# the node (Bioconductor's release may be older than the node needs), as the
# tag's source archive, which needs no GitHub API access or token; its
# dependencies come from Bioconductor and CRAN. lcmsPlot has no compiled code,
# so no build tools are needed. rawrr reads the Thermo .raw files; its
# rawrr.exe (which includes the Thermo RawFileReader) is downloaded once.
# An lcmsPlot already at that version or later is left as it is. The optional
# git ref (a branch or commit) replaces the tag, for testing a release before
# it is tagged.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || !nzchar(args[1])) {
    stop("Usage: Rscript install-packages.R <lcmsPlot version> [git ref]",
         call. = FALSE)
}
required <- args[1]
ref <- if (length(args) >= 2 && nzchar(args[2])) args[2] else paste0("v", required)
lib <- .libPaths()[1]

step <- function(...) cat("\n==> ", ..., "\n", sep = "")
fail <- function(...) {
    cat("\nERROR: ", ..., "\n", sep = "")
    quit(status = 1, save = "no")
}
installed_version <- function(pkg) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        as.character(utils::packageVersion(pkg))
    } else {
        NA_character_
    }
}

options(repos = c(CRAN = "https://cloud.r-project.org"))
cat("R ", as.character(getRversion()), "; library: ", lib, "\n", sep = "")

step("Installing BiocManager and remotes")
for (pkg in c("BiocManager", "remotes", "jsonlite")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        utils::install.packages(pkg, lib = lib)
    }
    if (!requireNamespace(pkg, quietly = TRUE)) fail("could not install ", pkg)
}
options(repos = BiocManager::repositories())

step("Installing lcmsPlot ", required, " (", ref, ")")
current <- installed_version("lcmsPlot")
if (!is.na(current) && utils::compareVersion(current, required) >= 0) {
    cat("lcmsPlot ", current, " is already installed.\n", sep = "")
} else {
    tryCatch(
        remotes::install_url(
            paste0("https://github.com/computational-metabolomics/lcmsPlot/",
                   "archive/", ref, ".tar.gz"),
            lib = lib, upgrade = "never", dependencies = NA,
            repos = BiocManager::repositories()),
        error = function(e) fail("could not install lcmsPlot: ",
                                 conditionMessage(e)))
    current <- installed_version("lcmsPlot")
    if (is.na(current) || utils::compareVersion(current, required) < 0) {
        fail("lcmsPlot ", required, " was not installed (found: ",
             if (is.na(current)) "none" else current, ").")
    }
}
cat("lcmsPlot ", current, "\n", sep = "")

step("Installing rawrr, which reads Thermo .raw files")
if (is.na(installed_version("rawrr"))) {
    BiocManager::install("rawrr", lib = lib, update = FALSE, ask = FALSE)
}
if (is.na(installed_version("rawrr"))) fail("could not install rawrr")
tryCatch(rawrr::installRawrrExe(),
         error = function(e) fail("could not download rawrr.exe: ",
                                  conditionMessage(e)))
if (!file.exists(rawrr:::.rawrrAssembly())) {
    fail("rawrr.exe was not found after the download.")
}

cat("\nAll R packages are installed.\n")
