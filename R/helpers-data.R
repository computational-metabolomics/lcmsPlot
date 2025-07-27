#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL

merge_by_index <- function(a, b, index_col) {
  b_mod <- b %>%
    mutate(row_id = row_number())

  a %>%
    left_join(b_mod, by = setNames("row_id", index_col))
}

remove_null_elements <- function(lst) {
  return(lst[!sapply(lst, is.null)])
}

get_mz_range <- function(mz, ppm = 5) {
  mzdev <- mz * (ppm / 1000000)
  return(c(mz - mzdev, mz + mzdev))
}

get_feature_data <- function(feature, options) {
  if (is.vector(feature)) {
    if (!is.na(feature['mz']) & !is.na(feature['rt'])) {
      feature_id <- paste0('M', round(feature['mz']), 'T', round(feature['rt']))
      mzr <- get_mz_range(feature['mz'], options$chromatograms$ppm)
      rtr <- c(feature['rt'] - options$chromatograms$rt_tol, feature['rt'] + options$chromatograms$rt_tol)
    } else if (!is.na(feature['mzmin']) & !is.na(feature['mzmax']) & !is.na(feature['rtmin']) & !is.na(feature['rtmax'])) {
      mzr <- c(feature['mzmin'], feature['mzmax'])
      rtr <- c(feature['rtmin'], feature['rtmax'])
      feature_id <- paste0('M', round((feature['mzmin'] + feature['mzmax']) / 2), 'T', round((feature['rtmin'] + feature['rtmax']) / 2))
    } else {
      stop("Features format not supported")
    }
  } else if (is.data.frame(feature)) {
    feature_id <- feature$name
    mzr <- get_mz_range(feature$mz, options$chromatograms$ppm)
    rtr <- c(feature$rt - options$chromatograms$rt_tol, feature$rt + options$chromatograms$rt_tol)
  } else {
    stop("Features format not supported")
  }

  return(list(
    feature_id = feature_id,
    mzr = mzr,
    rtr = rtr
  ))
}
