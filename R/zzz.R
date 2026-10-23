#' @importFrom methods show new validObject slot slotNames is as
#' @importFrom stats as.formula setNames median
#' @import dplyr
#' @import tibble
#' @import ggplot2
#' @importClassesFrom MsExperiment MsExperiment
#' @importClassesFrom Spectra Spectra
#' @importClassesFrom Chromatograms Chromatograms
NULL

utils::globalVariables(c(
    ".", ".data",
    # Columns referenced in dplyr::join_by() and data.table expressions.
    "sample_index", "mz", "mzMin", "mzMax", "rt", "rtMin", "rtMax",
    "intensity"
))

# bin_points() groups with data.table's `[` without importing data.table.
.datatable.aware <- TRUE
