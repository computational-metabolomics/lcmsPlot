#' Get the metadata associated with the input object
#'
#' `get_metadata()` is a generic helper used internally to standardise
#' metadata extraction across different classes of input objects.
#'
#' Depending on the class of `obj`, metadata may be:
#' - **constructed** (e.g., from character vectors of file paths), or
#' - **extracted and optionally replaced** (e.g., from `XCMSnExp` or
#'   `MsExperiment` objects).
#'
#' Across all methods, the returned metadata is enriched with:
#' - **`sample_index`** - a sequential index of samples
#' - **`sample_id`** - an identifier column selected via `sample_id_column`
#' - **`sample_path`** - a file path associated with each sample (if
#'   applicable)
#'
#' @section Character vector input (`character`):
#' A character vector is treated as a list of sample file paths (e.g.,
#' `.mzML`, `.mzXML`, `.cdf`).
#'
#' **If `metadata` is NULL:**
#' Metadata is *constructed automatically*:
#' * `sample_path`: full paths given in `obj`
#' * `sample_index`: row number
#' * `sample_id`: basename of each file without extension
#'
#' **If `metadata` is provided:**
#' The supplied `metadata` is used and the following columns are added:
#' * `sample_index`: row number
#' * `sample_id`: extracted from the `sample_id_column`
#' * `sample_path`: the input paths from `obj`
#'
#'
#' @section `XCMSnExp` input:
#' Metadata is taken from `xcms::phenoData(obj)`.
#'
#' **If `metadata` is provided:**
#' It replaces the existing `phenoData`.
#'
#' The returned metadata always includes:
#' * `sample_index`: row number
#' * `sample_id`: extracted using `sample_id_column`
#' * `sample_path`: values from `xcms::fileNames(obj)`
#'
#' @section `MsExperiment` input:
#' Metadata is taken from `MsExperiment::sampleData(obj)`.
#'
#' **If `metadata` is provided:**
#' It replaces existing `sampleData`.
#'
#' The returned metadata includes:
#' * `sample_index`: row number
#' * `sample_id`: extracted using `sample_id_column`
#' * `sample_path`: values from `xcms::fileNames(obj)`
#'
#' @section `MChromatograms` and `XChromatograms` input:
#' Metadata is taken from `MSnbase::phenoData(obj)`.
#'
#' `sample_id` is resolved in the following order:
#' 1. `sample_id_column` (if provided and present in `phenoData`)
#' 2. Basename without extension of the `spectraOrigin` column
#' 3. `"sample1"`, `"sample2"`, ... (fallback)
#'
#' `sample_path` is set to `spectraOrigin` if present, otherwise `NA`.
#'
#' **If `metadata` is provided:**
#' It is column-bound onto the extracted `phenoData`.
#'
#' The returned metadata always includes:
#' * `sample_index`: row number
#' * `sample_id`: resolved as above
#' * `sample_path`: from `spectraOrigin` or `NA`
#'
#' @section `XChromatogram` input:
#' Represents a single chromatogram (one sample). `sample_path` is always `NA`.
#'
#' **If `metadata` is NULL:**
#' Returns a single-row tibble with `sample_index = 1`, `sample_id = "sample1"`,
#' and `sample_path = NA`.
#'
#' **If `metadata` is provided:**
#' Must have exactly one row. The supplied metadata is used with:
#' * `sample_index`: `1`
#' * `sample_id`: from `sample_id_column` if present, otherwise `"sample1"`
#' * `sample_path`: `NA`
#'
#' @section `XcmsRawList` input:
#' Each element of the list is an `xcmsRaw` object. `sample_path` is extracted
#' from the `@@filepath` slot of each `xcmsRaw`.
#'
#' **If `metadata` is NULL:**
#' Metadata is *constructed automatically*:
#' * `sample_index`: sequential index
#' * `sample_id`: basename without extension of `sample_path`
#' * `sample_path`: from `xcmsRaw@@filepath`
#'
#' **If `metadata` is provided:**
#' Must have one row per `xcmsRaw` object. The supplied metadata is used with:
#' * `sample_index`: row number
#' * `sample_id`: from `sample_id_column` if present, otherwise basename of path
#' * `sample_path`: from `xcmsRaw@@filepath`
#'
#' @section `ExternalDataSource` input:
#' Metadata is taken from the `@@metadata` slot of the `ExternalDataSource`
#' object. The `metadata` parameter is ignored.
#'
#' The returned metadata always includes:
#' * `sample_index`: row number
#' * `sample_id`: extracted using `sample_id_column`
#'
#' @section `DBIConnection` input:
#' Metadata is queried from a Compound Discoverer SQLite database via
#' `get_workflow_input_files()`.
#'
#' **If `metadata` is NULL:**
#' Returns the query result with:
#' * `sample_index`: row number
#' * `sample_id`: from `StudyFileID`
#' * `sample_path`: from `PhysicalFileName`
#'
#' **If `metadata` is provided:**
#' It is joined onto the query result using `sample_id_column`.
#'
#' @param obj A data object containing or representing samples.
#' @param sample_id_column A `character` value indicating the column that
#' should be used as the sample ID.
#' @param metadata Optional metadata `tibble` used to replace or augment
#' sample metadata when not already embedded in the object.
#' @return A `tibble` containing standardised metadata with at least
#' `sample_index`, `sample_id`, and `sample_path`.
#' @keywords internal
get_metadata <- function(obj, sample_id_column, metadata) {
    UseMethod("get_metadata")
}

#' @rdname get_metadata
#' @keywords internal
#' @exportS3Method
get_metadata.character <- function(obj, sample_id_column, metadata) {
    if (is.null(metadata)) {
        tibble(sample_path = obj) |>
            mutate(
                sample_index = row_number(),
                sample_id = tools::file_path_sans_ext(basename(obj))
            )
    } else {
        as_tibble(metadata) |>
            mutate(
                sample_index = row_number(),
                sample_id = .data[[sample_id_column]],
                sample_path = obj
            )
    }
}

#' @rdname get_metadata
#' @keywords internal
#' @exportS3Method
get_metadata.XCMSnExp <- function(obj, sample_id_column, metadata) {
    if (!is.null(metadata)) {
        xcms::phenoData(obj) <- new("AnnotatedDataFrame", metadata)
    }

    paths <- xcms::fileNames(obj)
    df <- as_tibble(xcms::phenoData(obj)@data)
    df |>
        mutate(
            sample_index = row_number(),
            sample_id = if (
                !is.null(sample_id_column) && sample_id_column %in% colnames(df)
            ) {
                df[[sample_id_column]]
            } else {
                tools::file_path_sans_ext(basename(paths))
            },
            sample_path = paths
        )
}

#' @rdname get_metadata
#' @keywords internal
#' @exportS3Method
get_metadata.MsExperiment <- function(obj, sample_id_column, metadata) {
    if (!is.null(metadata)) {
        MsExperiment::sampleData(obj) <- metadata
    }

    paths <- xcms::fileNames(obj)
    df <- MsExperiment::sampleData(obj) |> as_tibble()
    df |>
        mutate(
            sample_index = row_number(),
            sample_id = if (
                !is.null(sample_id_column) && sample_id_column %in% colnames(df)
            ) {
                df[[sample_id_column]]
            } else {
                tools::file_path_sans_ext(basename(paths))
            },
            sample_path = paths
        )
}

#' @rdname get_metadata
#' @keywords internal
#' @exportS3Method
get_metadata.MChromatograms <- function(obj, sample_id_column, metadata) {
    df <- MSnbase::phenoData(obj)@data |> as_tibble()

    if (!is.null(sample_id_column) && sample_id_column %in% colnames(df)) {
        sample_ids <- df[[sample_id_column]]
    } else if ("spectraOrigin" %in% colnames(df)) {
        sample_ids <- tools::file_path_sans_ext(basename(df$spectraOrigin))
    } else {
        sample_ids <- paste0("sample", seq_len(nrow(df)))
    }

    sample_paths <- if ("spectraOrigin" %in% colnames(df)) {
        df$spectraOrigin
    } else {
        rep(NA_character_, nrow(df))
    }

    df_metadata <- df |>
        mutate(
            sample_index = row_number(),
            sample_id = sample_ids,
            sample_path = sample_paths
        )

    if (!is.null(metadata)) {
        df_metadata <- as_tibble(cbind(df_metadata, metadata))
    }

    return(df_metadata)
}

#' @rdname get_metadata
#' @keywords internal
#' @exportS3Method
get_metadata.XChromatograms <- function(obj, sample_id_column, metadata) {
    df <- MSnbase::phenoData(obj)@data |> as_tibble()

    if (!is.null(sample_id_column) && sample_id_column %in% colnames(df)) {
        sample_ids <- df[[sample_id_column]]
    } else if ("spectraOrigin" %in% colnames(df)) {
        sample_ids <- tools::file_path_sans_ext(basename(df$spectraOrigin))
    } else {
        sample_ids <- paste0("sample", seq_len(nrow(df)))
    }

    sample_paths <- if ("spectraOrigin" %in% colnames(df)) {
        df$spectraOrigin
    } else {
        rep(NA_character_, nrow(df))
    }

    df_metadata <- df |>
        mutate(
            sample_index = row_number(),
            sample_id = sample_ids,
            sample_path = sample_paths
        )

    if (!is.null(metadata)) {
        df_metadata <- as_tibble(cbind(df_metadata, metadata))
    }

    return(df_metadata)
}

#' @rdname get_metadata
#' @keywords internal
#' @exportS3Method
get_metadata.XChromatogram <- function(obj, sample_id_column, metadata) {
    if (!is.null(metadata)) {
        if (nrow(metadata) != 1) {
            stop("metadata needs to have one element for XChromatogram input.")
        }

        as_tibble(metadata) |>
            mutate(
                sample_index = row_number(),
                sample_id = if (
                    !is.null(sample_id_column) &&
                    sample_id_column %in% colnames(metadata)
                ) {
                    .data[[sample_id_column]]
                } else {
                    paste0("sample", row_number())
                },
                sample_path = NA_character_
            )
    } else {
        tibble(
            sample_index = 1L,
            sample_id = "sample1",
            sample_path = NA_character_
        )
    }
}

#' @rdname get_metadata
#' @keywords internal
#' @exportS3Method
get_metadata.XcmsRawList <- function(obj, sample_id_column, metadata) {
    objs <- obj@data
    paths <- vapply(objs, function(x) {
        if (length(x@filepath) > 0L) x@filepath[[1L]] else NA_character_
    }, character(1L))
    default_ids <- tools::file_path_sans_ext(basename(paths))

    if (!is.null(metadata)) {
        if (nrow(metadata) != length(objs)) {
            stop(
                "metadata must have one row per xcmsRaw object (",
                length(objs), " expected, ", nrow(metadata), " provided)."
            )
        }

        as_tibble(metadata) |>
            mutate(
                sample_index = row_number(),
                sample_id = if (
                    !is.null(sample_id_column) &&
                    sample_id_column %in% colnames(metadata)
                ) {
                    .data[[sample_id_column]]
                } else {
                    default_ids
                },
                sample_path = paths
            )
    } else {
        tibble(
            sample_index = seq_along(objs),
            sample_id = default_ids,
            sample_path = paths
        )
    }
}

#' @rdname get_metadata
#' @keywords internal
#' @exportS3Method
get_metadata.ExternalDataSource <- function(obj, sample_id_column, metadata) {
    obj@metadata |>
        as_tibble() |>
        mutate(
            sample_index = row_number(),
            sample_id = .data[[sample_id_column]]
        )
}

#' @rdname get_metadata
#' @keywords internal
#' @exportS3Method
get_metadata.DBIConnection <- function(obj, sample_id_column, metadata) {
    cd_metadata <- get_workflow_input_files(obj) |>
        dplyr::mutate(
            sample_index = dplyr::row_number(),
            sample_id = .data$StudyFileID,
            sample_path = .data$PhysicalFileName
        )

    if (is.null(metadata)) {
        return(cd_metadata)
    }

    if (!sample_id_column %in% colnames(metadata)) {
        stop(sprintf(
            "Column '%s' not found in metadata",
            sample_id_column
        ))
    }

    dplyr::left_join(
        metadata,
        cd_metadata,
        by = setNames("sample_id", sample_id_column)
    )
}

#' @rdname get_metadata
#' @keywords internal
#' @exportS3Method
get_metadata.purityA <- function(obj, sample_id_column, metadata) {
    paths <- obj@fileList
    default_ids <- tools::file_path_sans_ext(basename(paths))

    if (!is.null(metadata)) {
        as_tibble(metadata) |>
            mutate(
                sample_index = row_number(),
                sample_id = if (
                    !is.null(sample_id_column) &&
                    sample_id_column %in% colnames(metadata)
                ) {
                    .data[[sample_id_column]]
                } else {
                    default_ids
                },
                sample_path = paths
            )
    } else {
        tibble(
            sample_index = seq_along(paths),
            sample_id = default_ids,
            sample_path = paths
        )
    }
}

#' Get the detected peaks from the data object (e.g. XCMSnExp)
#'
#' `get_detected_peaks()` is an internal helper that standardises extraction of
#' detected chromatographic peaks across different object types commonly used in
#' LC-MS workflows.
#'
#' Supported inputs behave as follows:
#'
#' - **`character`** – Assumed to represent sample paths; no peak detection
#' information is available. Always returns `NULL`.
#'
#' - **`XCMSnExp`** and **`MsExperiment`** – If the object is processed and
#' contains chromatographic peaks, extracts `xcms::chromPeaks(obj)` and
#' returns it as a data frame.
#' The column `sample` is renamed to `sample_index`.
#'
#' When peaks are not found or the object is not processed, `NULL` is returned.
#'
#' @param obj A data object containing or representing samples.
#' @return A `tibble` of detected peaks (one row per peak), or `NULL` if no
#' peaks are available.
#' @keywords internal
get_detected_peaks <- function(obj) {
    UseMethod("get_detected_peaks")
}

.get_detected_peaks_xcms <- function(obj) {
    if (is_xcms_processed_data(obj) && xcms::hasChromPeaks(obj)) {
        as_tibble(xcms::chromPeaks(obj)) |>
            dplyr::rename(sample_index = sample)
    } else {
        NULL
    }
}

#' @rdname get_detected_peaks
#' @keywords internal
#' @exportS3Method
get_detected_peaks.character <- function(obj) {
    NULL
}

#' @rdname get_detected_peaks
#' @keywords internal
#' @exportS3Method
get_detected_peaks.XCMSnExp <- function(obj) {
    .get_detected_peaks_xcms(obj)
}

#' @rdname get_detected_peaks
#' @keywords internal
#' @exportS3Method
get_detected_peaks.MsExperiment <- function(obj) {
    .get_detected_peaks_xcms(obj)
}

#' @rdname get_detected_peaks
#' @keywords internal
#' @exportS3Method
get_detected_peaks.MChromatograms <- function(obj) {
    NULL
}

#' @rdname get_detected_peaks
#' @keywords internal
#' @exportS3Method
get_detected_peaks.XChromatograms <- function(obj) {
    if (any(xcms::hasChromPeaks(obj))) {
        peaks <- as_tibble(xcms::chromPeaks(obj)) |>
            dplyr::rename(sample_index = "column") |>
            select(-dplyr::all_of("mz")) # This is already present in mz_info

        # Extract mz ranges from each row of the XChromatograms object
        mz_info <- do.call(rbind, lapply(seq_len(nrow(obj)), function(i) {
            mz_range <- MSnbase::mz(obj[i, 1L])[[1L]]
            tibble(
                row   = i,
                mz    = mean(mz_range)
            )
        }))

        dplyr::left_join(peaks, mz_info, by = "row")
    } else {
        NULL
    }
}

#' @rdname get_detected_peaks
#' @keywords internal
#' @exportS3Method
get_detected_peaks.XChromatogram <- function(obj) {
    if (xcms::hasChromPeaks(obj)) {
        mz_range <- MSnbase::mz(obj)
        as_tibble(xcms::chromPeaks(obj)) |>
            mutate(
                mz = mean(mz_range),
                sample_index = obj@fromFile
            )
    } else {
        NULL
    }
}

#' @rdname get_detected_peaks
#' @keywords internal
#' @exportS3Method
get_detected_peaks.XcmsRawList <- function(obj) {
    NULL
}

#' @rdname get_detected_peaks
#' @keywords internal
#' @exportS3Method
get_detected_peaks.ExternalDataSource <- function(obj) {
    obj@peaks
}

#' @rdname get_detected_peaks
#' @keywords internal
#' @exportS3Method
get_detected_peaks.purityA <- function(obj) {
    if (is.data.frame(obj@grped_df) && nrow(obj@grped_df) > 0) {
        as_tibble(obj@grped_df) |>
            dplyr::rename(sample_index = "sample") |>
            dplyr::select(
                dplyr::any_of(
                    c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax",
                      "into", "maxo", "sample_index")
                )
            )
    } else {
        NULL
    }
}

#' Get the grouped peaks across samples (features) from the data object
#'
#' `get_grouped_peaks()` is an internal helper that retrieves feature-level
#' grouped peaks, i.e., chromatographic peaks aligned across samples.
#'
#' @param obj A data object containing or representing samples.
#' @return A `tibble` of grouped (feature-level) peaks,
#' or `NULL` if not available.
#' @keywords internal
get_grouped_peaks <- function(obj) {
    UseMethod("get_grouped_peaks")
}

#' @rdname get_grouped_peaks
#' @keywords internal
#' @exportS3Method
get_grouped_peaks.default <- function(obj) {
    return(NULL)
}

.get_grouped_peaks_xcms <- function(obj) {
    if (is_xcms_processed_data(obj) && xcms::hasFeatures(obj)) {
        as_tibble(xcms::featureDefinitions(obj)) |>
            rename(all_of(c(mz = "mzmed", rt = "rtmed"))) |>
            mutate(name = xcms_utils$group_names(obj)) |>
            xcms_utils$format_feature_identifiers(
                num_digits_rt = 0,
                num_digits_mz = 4)
    } else {
        NULL
    }
}

#' @rdname get_grouped_peaks
#' @keywords internal
#' @exportS3Method
get_grouped_peaks.XCMSnExp <- function(obj) {
    .get_grouped_peaks_xcms(obj)
}

#' @rdname get_grouped_peaks
#' @keywords internal
#' @exportS3Method
get_grouped_peaks.XcmsExperiment <- function(obj) {
    .get_grouped_peaks_xcms(obj)
}

#' @rdname get_grouped_peaks
#' @keywords internal
#' @exportS3Method
get_grouped_peaks.MsExperiment <- function(obj) {
    .get_grouped_peaks_xcms(obj)
}

#' @rdname get_grouped_peaks
#' @keywords internal
#' @exportS3Method
get_grouped_peaks.purityA <- function(obj) {
    if (is.data.frame(obj@grped_df) && nrow(obj@grped_df) > 0) {
        as_tibble(obj@grped_df) |>
            dplyr::group_by(.data$grpid) |>
            dplyr::summarise(
                mz = mean(.data$mz, na.rm = TRUE),
                rt = mean(.data$rt, na.rm = TRUE),
                .groups = "drop"
            ) |>
            dplyr::mutate(name = as.character(.data$grpid))
    } else {
        NULL
    }
}
