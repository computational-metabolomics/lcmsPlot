#' Construct a field specification for data frame validation
#'
#' @param name A `character` value giving the column name.
#' @param type A function used to check the column's type (e.g. `is.numeric`).
#' If `NULL`, no type check is performed.
#' @param required A `logical` value indicating whether the column must be
#' present. Defaults to `TRUE`.
#' @return A named `list` with elements `name`, `type`, and `required`.
#' @keywords internal
field <- function(name, type = NULL, required = TRUE) {
    list(name = name, type = type, required = required)
}

#' Validate a data frame against a list of field specifications
#'
#' Returns `TRUE` on success or a `character` string describing the first
#' failure. `NULL` and empty data frames always pass.
#'
#' @param df A `data.frame` or `tibble` to validate.
#' @param fields A `list` of field specifications created with `field()`.
#' @param exact A `logical` value. If `TRUE`, the column names of `df` must
#' be identical to the names in `fields` (same set, same order). Defaults
#' to `FALSE`.
#' @return `TRUE` if validation passes, or a `character` string describing
#' the first failure.
#' @keywords internal
validate_data_frame <- function(df, fields, exact = FALSE) {
    if (is.null(df) || nrow(df) == 0) {
        return(TRUE)
    }

    if (!is_tibble(df)) {
        return("The data is not a tibble.")
    }

    for (field in fields) {
        present <- field$name %in% colnames(df)
        if (field$required && !present) {
            return(sprintf("missing required column '%s'", field$name))
        }
        if (present && !is.null(field$type) && !field$type(df[[field$name]])) {
            return(sprintf("column '%s' has unexpected type", field$name))
        }
    }

    if (exact) {
        expected <- vapply(fields, `[[`, "name", FUN.VALUE = character(1))
        if (!identical(colnames(df), expected)) {
            return(sprintf(
                "expected columns [%s], got [%s]",
                paste(expected, collapse = ", "),
                paste(colnames(df), collapse = ", ")
            ))
        }
    }

    TRUE
}

#' Validate the data frame slots of an S4 object
#'
#' Iterates over a named list of validator functions, applying each to the
#' corresponding slot of `object`.
#'
#' @param object An S4 object whose slots are to be validated.
#' @param validators A named `list` of functions, where each name matches a
#' slot of `object` and each function accepts a `data.frame`/`tibble` and returns
#' `TRUE` or a character error string (as produced by `validate_data_frame()`).
#' @return `TRUE` if all validators pass, or a `character` string describing
#' the first failure in the form `"@slot_name: <message>"`.
#' @keywords internal
validate_object <- function(object, validators) {
    ret <- TRUE

    for (validator_name in names(validators)) {
        df <- slot(object, validator_name)
        result <- validators[[validator_name]](df)

        if (!isTRUE(result)) {
            ret <- sprintf("@%s: %s", validator_name, result)
            break
        }
    }

    ret
}
