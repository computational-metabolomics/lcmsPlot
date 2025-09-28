#' Get the mzR objects for a set of sample paths.
#'
#' @param sample_paths The sample paths.
#' @returns A list of mzR objects.
io_get_raw_data <- function(sample_paths) {
  raw_data <- list()

  for (sample_path in sample_paths) {
    raw_data[[sample_path]] <- mzR::openMSfile(sample_path)
  }

  return(raw_data)
}

#' Close a list of mzR object connections.
#'
#' @param raw_data A list of mzR object connections.
io_close_raw_data <- function(raw_data) {
  for (raw_obj in raw_data) {
    mzR::close(raw_obj)
  }
}
