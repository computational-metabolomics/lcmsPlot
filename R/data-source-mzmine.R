#' Create an MZmine data source from feature lists
#'
#' This function reads MZmine feature list files (supports version 2 and above)
#' and sample metadata and converts them into a format compatible
#' with xcms-style peak tables.
#'
#' A feature list may hold one data file (a list exported per sample) or
#' several (a list aligned across samples, e.g. after the join aligner). In an
#' aligned list, each data file contributes the peaks it has: rows whose state
#' in that file is `UNKNOWN`, or which have no retention time there, are
#' skipped, while gap-filled (`ESTIMATED`) peaks are kept.
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
#' mzmine_features <- tibble::tibble(
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

        # One feature state/status column per data file in the list.
        state_cols <- if (is_mzmine2) {
            grep(" Feature status$", colnames(df), value = TRUE)
        } else {
            grep("^datafile:.*?:feature_state$", colnames(df), value = TRUE)
        }

        if (length(state_cols) == 0) {
            stop("Expected a feature state/status column for each data file")
        }

        filenames <- if (is_mzmine2) {
            sub(" Feature status$", "", state_cols)
        } else {
            sub("^datafile:(.*):feature_state$", "\\1", state_cols)
        }

        do.call(rbind, Map(function(state_col, filename) {
            .mzmine_file_peaks(
                df, state_col, filename, is_mzmine2,
                single_file = length(state_cols) == 1,
                sample_basenames = sample_basenames)
        }, state_cols, filenames, USE.NAMES = FALSE))
    })

    peaks <- do.call(rbind, peaks_list) |>
        as_tibble() |>
        convert_rt_to_seconds()

    new(
        "ExternalDataSource",
        name = "mzmine",
        metadata = metadata,
        peaks = peaks
    )
}

#' The peaks of one data file of an MZmine feature list
#'
#' @param df A `data.frame` holding the feature list.
#' @param state_col The name of the data file's feature state/status column.
#' @param filename The data file's name, as used in the column names.
#' @param is_mzmine2 A `logical` value indicating an MZmine 2 feature list.
#' @param single_file A `logical` value indicating that the list holds only
#' this data file; an MZmine 3+ list exported per sample then keeps the peak
#' values in its row-level columns.
#' @param sample_basenames A `character` vector of the sample file names.
#' @return A `data.frame` of the file's peaks, with the xcms-style columns
#' `mz`, `rt`, `into`, `maxo`, `mzmin`, `mzmax`, `rtmin`, `rtmax` and
#' `sample_index`.
#' @keywords internal
.mzmine_file_peaks <- function(
    df,
    state_col,
    filename,
    is_mzmine2,
    single_file,
    sample_basenames
) {
    sample_index <- match(basename(filename), sample_basenames)
    if (is.na(sample_index)) {
        stop("Could not match feature list file to sample_paths")
    }

    targets <- c("mz", "rt", "into", "maxo", "mzmin", "mzmax", "rtmin", "rtmax")

    if (is_mzmine2) {
        col_map <- setNames(targets, paste0(filename, c(
            " Feature m/z", " Feature RT", " Peak area", " Peak height",
            " Feature m/z min", " Feature m/z max",
            " Feature RT start", " Feature RT end")))
    } else {
        fields <- c("mz", "rt", "area", "height", "mz_range:min",
                    "mz_range:max", "rt_range:min", "rt_range:max")
        per_file <- paste0("datafile:", filename, ":", fields)
        col_map <- if (all(per_file %in% colnames(df)) || !single_file) {
            setNames(targets, per_file)
        } else {
            setNames(targets, fields)
        }
    }

    missing <- setdiff(names(col_map), colnames(df))
    if (length(missing)) {
        stop("Missing required columns in feature list")
    }

    out <- df[, names(col_map), drop = FALSE]
    colnames(out) <- col_map

    # Skip the rows this data file has no peak for.
    state <- as.character(df[[state_col]])
    has_peak <- !is.na(out$rt) & !(state %in% c("UNKNOWN", "") | is.na(state))
    out <- out[has_peak, , drop = FALSE]

    out$sample_index <- rep(sample_index, nrow(out))
    out
}
