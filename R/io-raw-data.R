#' Get the `mzR` objects for a set of sample paths
#'
#' @param sample_paths A `character` vector of file paths to the
#' raw MS data files (e.g., `.mzML`, `.mzXML`, `.CDF`).
#' @return A named list of `mzR` objects, where each element corresponds to a
#' file in `sample_paths`.
#' @keywords internal
io_get_raw_data <- function(sample_paths) {
    raw_data <- list()

    for (sample_path in sample_paths) {
        raw_data[[sample_path]] <- mzR::openMSfile(sample_path)
    }

    return(raw_data)
}

#' Close open `mzR` object connections
#'
#' @param raw_data A list of `mzR` objects (as returned by `io_get_raw_data()`).
#' @return `NULL`
#' @keywords internal
io_close_raw_data <- function(raw_data) {
    for (raw_obj in raw_data) {
        mzR::close(raw_obj)
    }
}
