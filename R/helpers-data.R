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

#' @import dplyr
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
