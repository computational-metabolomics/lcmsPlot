#' Create an MS-DIAL data source from peak lists
#'
#' This function reads MS-DIAL peak list files and sample metadata
#' and converts them into a format compatible with xcms-style peak tables.
#'
#' @param peaks_paths A `character` vector of paths to MS-DIAL
#' peak list files (CSV or TSV). The ordering should match the one
#' in the metadata (i.e., the samples).
#' @param sample_paths A `character` vector of paths to raw sample files
#' (e.g., `mzML`).
#' @param metadata_path A `character` value (optional) indicating the
#' path to a metadata file (CSV or TSV).
#' If `NULL`, a metadata `data.frame` is created from `sample_paths`.
#' @return An object of class `ExternalDataSource`.
#' @export
MsDialPeaksSource <- function(
    peaks_paths,
    sample_paths,
    metadata_path = NULL
) {
    metadata <- process_metadata(sample_paths, metadata_path)

    peaks_list <- lapply(seq_along(peaks_paths), function(i) {
        path <- peaks_paths[i]
        df <- utils::read.table(
            path,
            sep = detect_separator(path),
            header = TRUE,
            stringsAsFactors = FALSE,
            check.names = FALSE
        )

        sample_index <- i

        if (sample_index > length(sample_paths)) {
            stop("Number of peak lists exceeds number of sample paths")
        }

        # Column mapping
        col_map <- c(
            "Precursor m/z" = "mz",
            "RT (min)" = "rt",
            "Area" = "into",
            "Height" = "maxo",
            "RT left(min)" = "rtmin",
            "RT right (min)" = "rtmax"
        )

        missing <- setdiff(names(col_map), colnames(df))
        if (length(missing)) {
            stop(paste(
                "Missing required columns in MS-DIAL peak list:",
                paste(missing, collapse = ", ")
            ))
        }

        out <- df[, names(col_map), drop = FALSE]
        colnames(out) <- col_map

        out$sample_index <- sample_index
        out
    })

    peaks <- do.call(rbind, peaks_list) |>
        convert_rt_to_seconds()

    new(
        "ExternalDataSource",
        name = "ms-dial",
        metadata = metadata,
        peaks = peaks
    )
}
