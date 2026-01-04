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
