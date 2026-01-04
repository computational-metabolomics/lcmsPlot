#' Create an MZmine data source from feature lists
#'
#' This function reads MZmine feature list files (supports version 2 and above)
#' and sample metadata and converts them into a format compatible
#' with xcms-style peak tables.
#'
#' @param feature_lists_paths A `character` vector of paths to MZmine
#' feature list files (CSV or TSV).
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
#' ## Fake sample paths
#' sample_paths <- file.path(
#'   tmp_dir,
#'   c("sample1.mzML", "sample2.mzML")
#' )
#'
#' ## Create a minimal MZmine 2 feature list CSV
#' feature_list_path <- file.path(tmp_dir, "mzmine_features.csv")
#'
#' mzmine_features <- data.frame(
#'   "row m/z" = c(100.1, 200.2),
#'   "sample1.mzML Feature status" = c("DETECTED", "DETECTED"),
#'   "sample1.mzML Feature m/z" = c(100.1, 200.2),
#'   "sample1.mzML Feature RT" = c(300, 600),
#'   "sample1.mzML Peak area" = c(10000, 20000),
#'   "sample1.mzML Peak height" = c(500, 800),
#'   "sample1.mzML Feature m/z min" = c(99.9, 199.9),
#'   "sample1.mzML Feature m/z max" = c(100.3, 200.4),
#'   "sample1.mzML Feature RT start" = c(290, 590),
#'   "sample1.mzML Feature RT end" = c(310, 610),
#'   check.names = FALSE
#' )
#'
#' utils::write.csv(mzmine_features, feature_list_path, row.names = FALSE)
#'
#' ## Create the data source
#' ds <- MZmineFeatureListsSource(
#'   feature_lists_paths = feature_list_path,
#'   sample_paths = sample_paths
#' )
MZmineFeatureListsSource <- function(
    feature_lists_paths,
    sample_paths,
    metadata_path = NULL
) {
    # 1. Metadata processing
    metadata <- process_metadata(sample_paths, metadata_path)
    sample_basenames <- basename(sample_paths)

    # 2. Feature lists processing
    peaks_list <- lapply(feature_lists_paths, function(path) {
        df <- utils::read.table(
            path,
            sep = detect_separator(path),
            header = TRUE,
            stringsAsFactors = FALSE,
            check.names = FALSE
        )

        is_mzmine2 <- "row m/z" %in% colnames(df)

        if (is_mzmine2) {
            state_col <- grep(" Feature status$", colnames(df), value = TRUE)
        } else {
            state_col <- grep("^datafile:.*?:feature_state$", colnames(df), value = TRUE)
        }

        if (length(state_col) != 1) {
            stop("Expected exactly one feature state/status column")
        }

        # Filename extraction
        filename <- if (is_mzmine2) {
            sub(" Feature status$", "", state_col)
        } else {
            sub("^datafile:(.*):feature_state$", "\\1", state_col)
        }

        filename <- basename(filename)

        sample_index <- match(filename, sample_basenames)
        if (is.na(sample_index)) {
            stop("Could not match feature list file to sample_paths")
        }

        # Column mapping
        if (is_mzmine2) {
            col_map <- setNames(
                c(
                    "mz", "rt", "into", "maxo",
                    "mzmin", "mzmax", "rtmin", "rtmax"
                ),
                c(
                    paste0(filename, " Feature m/z"),
                    paste0(filename, " Feature RT"),
                    paste0(filename, " Peak area"),
                    paste0(filename, " Peak height"),
                    paste0(filename, " Feature m/z min"),
                    paste0(filename, " Feature m/z max"),
                    paste0(filename, " Feature RT start"),
                    paste0(filename, " Feature RT end")
                )
            )
        } else {
            col_map <- c(
                "mz" = "mz",
                "rt" = "rt",
                "area" = "into",
                "height" = "maxo",
                "mz_range:min" = "mzmin",
                "mz_range:max" = "mzmax",
                "rt_range:min" = "rtmin",
                "rt_range:max" = "rtmax"
            )
        }

        missing <- setdiff(names(col_map), colnames(df))
        if (length(missing)) {
            stop("Missing required columns in feature list")
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
        name = "mzmine",
        metadata = metadata,
        peaks = peaks
    )
}
