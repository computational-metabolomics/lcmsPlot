#' External data source wrapper
#'
#' An S4 class representing an external data source such as MZmine.
#' This class acts as an adapter for external data sources by converting
#' data to a common format (XCMS).
#'
#' @slot name A `character` value indicating the name of the data source.
#' @slot metadata A `data.frame` representing the sample metadata.
#' @slot peaks A `data.frame` representing the exported peaks.
#' @export
setClass(
    "ExternalDataSource",
    slots = list(
        name = "character",
        metadata = "data.frame",
        peaks = "data.frame"
    ),
    prototype = list(
        name = NULL,
        metadata = NULL,
        peaks = NULL
    )
)

setValidity("ExternalDataSource", function(object) {
    if (is.null(object@name) | !is.character(object@name)) {
        return("name must be a string")
    }

    if (!is.data.frame(object@metadata)) {
        return("metadata must be a data.frame")
    }

    if (!"sample_path" %in% colnames(object@metadata)) {
        return("metadata must contain a 'sample_path' column")
    }

    if (is.null(object@peaks)) {
        return(TRUE)
    }

    if (!is.data.frame(object@peaks)) {
        return("peaks must be a data.frame")
    }

    required_cols <- c(
        "mz", "rt", "rtmin", "rtmax", "into", "maxo", "sample_index"
    )

    missing_cols <- setdiff(required_cols, colnames(object@peaks))
    if (length(missing_cols) > 0) {
        return(
            paste(
                "peaks is missing required columns:",
                paste(missing_cols, collapse = ", ")
            )
        )
    }

    TRUE
})

#' Show a summary of an instance of class `ExternalDataSource`
#'
#' @param object An instance of class `ExternalDataSource`.
#' @return Invisible \code{NULL}
#' @export
#' @examples
#' ## Create dummy metadata
#' metadata <- data.frame(
#'   sample_id = c("S1", "S2"),
#'   sample_path = c("sample1.mzML", "sample2.mzML"),
#'   stringsAsFactors = FALSE
#' )
#'
#' ## Create dummy peaks
#' peaks <- data.frame(
#'   mz = c(100.1, 150.2),
#'   rt = c(300, 450),
#'   rtmin = c(290, 440),
#'   rtmax = c(310, 460),
#'   into = c(10000, 15000),
#'   maxo = c(2000, 2500),
#'   sample_index = c(1, 2)
#' )
#'
#' ## Create ExternalDataSource object
#' eds <- new(
#'   "ExternalDataSource",
#'   name = "Example data source",
#'   metadata = metadata,
#'   peaks = peaks
#' )
#'
#' ## Show summary
#' eds
setMethod(
    f = "show",
    signature = "ExternalDataSource",
    function(object) {
        cat("Object of class", class(object), "\n")

        # Name
        cat(" Name:", object@name, "\n")

        # Metadata summary
        if (!is.null(object@metadata) && nrow(object@metadata) > 0) {
            cat(
                " Metadata:",
                paste0(
                    nrow(object@metadata),
                    " rows, ",
                    ncol(object@metadata),
                    " columns"
                ),
                "\n"
            )

            if ("sample_path" %in% colnames(object@metadata)) {
                cat(" Sample paths (first 3):\n")
                cat(
                    paste0("  ", utils::head(object@metadata$sample_path, 3), collapse = "\n"),
                    "\n"
                )
            }
        } else {
            cat(" Metadata: <empty>\n")
        }

        # Peaks summary
        if (!is.null(object@peaks) && nrow(object@peaks) > 0) {
            cat(
                " Peaks:",
                paste0(
                    nrow(object@peaks),
                    " rows, ",
                    ncol(object@peaks),
                    " columns"
                ),
                "\n"
            )
        } else {
            cat(" Peaks: <none>\n")
        }
    }
)

#' Process sample metadata
#'
#' Construct a metadata data frame aligned with a set of sample paths. If a
#' metadata file is provided, it is read from disk and combined with the sample
#' paths. Otherwise, a minimal metadata table is created.
#'
#' @param sample_paths A `character` vector of sample file paths.
#' @param metadata_path A `character` value indicating the path to a delimited
#' metadata file with a header.
#' @return A `data.frame` containing a `sample_path` column and
#' any additional metadata.
#' @keywords internal
process_metadata <- function(sample_paths, metadata_path = NULL) {
    if (!is.null(metadata_path)) {
        metadata <- utils::read.table(
            metadata_path,
            sep = detect_separator(metadata_path),
            header = TRUE,
            stringsAsFactors = FALSE
        )

        if (nrow(metadata) != length(sample_paths)) {
            stop("Metadata rows must match length of sample_paths")
        }

        metadata$sample_path <- sample_paths
    } else {
        metadata <- data.frame(
            sample_path = sample_paths,
            stringsAsFactors = FALSE
        )
    }

    metadata
}

#' Convert retention times to seconds
#'
#' Convert retention time columns from minutes to seconds.
#'
#' @param peaks A `data.frame` containing `rt`, `rtmin`, and `rtmax`
#' columns in minutes.
#'
#' @return `peaks` with retention time columns converted to seconds.
#' @keywords internal
convert_rt_to_seconds <- function(peaks) {
    peaks |> mutate(
        rt = .data$rt * 60,
        rtmin = .data$rtmin * 60,
        rtmax = .data$rtmax * 60
    )
}
