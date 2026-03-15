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
#' @examples
#' ## Create temporary example files
#' tmp_dir <- tempdir()
#'
#' ## Fake raw sample paths
#' sample_paths <- file.path(
#'   tmp_dir,
#'   c("sample1.mzML", "sample2.mzML")
#' )
#'
#' ## Create minimal MS-DIAL peak list files (one per sample)
#' peak_list_paths <- file.path(
#'   tmp_dir,
#'   c("sample1_peaks.csv", "sample2_peaks.csv")
#' )
#'
#' msdial_peaks_1 <- tibble::tibble(
#'   "Precursor m/z" = c(100.1, 200.2),
#'   "RT (min)" = c(5.0, 10.0),
#'   "Area" = c(10000, 20000),
#'   "Height" = c(500, 800),
#'   "RT left(min)" = c(4.8, 9.8),
#'   "RT right (min)" = c(5.2, 10.2)
#' )
#'
#' msdial_peaks_2 <- tibble::tibble(
#'   "Precursor m/z" = c(150.3, 250.4),
#'   "RT (min)" = c(6.0, 12.0),
#'   "Area" = c(15000, 25000),
#'   "Height" = c(600, 900),
#'   "RT left(min)" = c(5.8, 11.8),
#'   "RT right (min)" = c(6.2, 12.2)
#' )
#'
#' utils::write.csv(msdial_peaks_1, peak_list_paths[1], row.names = FALSE)
#' utils::write.csv(msdial_peaks_2, peak_list_paths[2], row.names = FALSE)
#'
#' ## Create the data source
#' ds <- MsDialPeaksSource(
#'   peaks_paths = peak_list_paths,
#'   sample_paths = sample_paths
#' )
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
        as_tibble() |>
        convert_rt_to_seconds()

    new(
        "ExternalDataSource",
        name = "ms-dial",
        metadata = metadata,
        peaks = peaks
    )
}
