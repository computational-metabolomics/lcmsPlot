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
#' @param obj A data object containing or representing samples.
#' @param sample_id_column A `character` value indicating the column that
#' should be used as the sample ID.
#' @param metadata Optional metadata `data.frame` used to replace or augment
#' sample metadata when not already embedded in the object.
#' @return A `data.frame` containing standardised metadata with at least
#' `sample_index`, `sample_id`, and `sample_path`.
#' @keywords internal
get_metadata <- function(obj, sample_id_column, metadata) {
    UseMethod("get_metadata")
}

#' @rdname get_metadata
#' @keywords internal
get_metadata.character <- function(obj, sample_id_column, metadata) {
    if (is.null(metadata)) {
        data.frame(sample_path = obj) |>
            mutate(
                sample_index = row_number(),
                sample_id = tools::file_path_sans_ext(basename(obj))
            )
    } else {
        metadata |>
            mutate(
                sample_index = row_number(),
                sample_id = .data[[sample_id_column]],
                sample_path = obj
            )
    }
}

#' @rdname get_metadata
#' @keywords internal
get_metadata.XCMSnExp <- function(obj, sample_id_column, metadata) {
    if (!is.null(metadata)) {
        xcms::phenoData(obj) <- new("AnnotatedDataFrame", metadata)
    }

    paths <- xcms::fileNames(obj)
    df <- xcms::phenoData(obj)@data
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
get_metadata.MsExperiment <- function(obj, sample_id_column, metadata) {
    if (!is.null(metadata)) {
        MsExperiment::sampleData(obj) <- metadata
    }

    paths <- xcms::fileNames(obj)
    df <- MsExperiment::sampleData(obj) |> as.data.frame()
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
get_metadata.MChromatograms <- function(obj, sample_id_column, metadata) {
    df <- MSnbase::phenoData(obj)@data |> as.data.frame()

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
        df_metadata <- cbind(df_metadata, metadata)
    }

    return(df_metadata)
}

#' @rdname get_metadata
#' @keywords internal
get_metadata.XChromatograms <- function(obj, sample_id_column, metadata) {
    df <- MSnbase::phenoData(obj)@data |> as.data.frame()

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
        df_metadata <- cbind(df_metadata, metadata)
    }

    return(df_metadata)
}

#' @rdname get_metadata
#' @keywords internal
get_metadata.XChromatogram <- function(obj, sample_id_column, metadata) {
    if (!is.null(metadata)) {
        if (nrow(metadata) != 1) {
            stop("metadata needs to have one element for XChromatogram input.")
        }

        metadata |>
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
        data.frame(
            sample_index = 1L,
            sample_id = "sample1",
            sample_path = NA_character_
        )
    }
}

#' @rdname get_metadata
#' @keywords internal
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

        metadata |>
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
        data.frame(
            sample_index = seq_along(objs),
            sample_id = default_ids,
            sample_path = paths
        )
    }
}

#' @rdname get_metadata
#' @keywords internal
get_metadata.ExternalDataSource <- function(obj, sample_id_column, metadata) {
    obj@metadata |>
        as.data.frame() |>
        mutate(
            sample_index = row_number(),
            sample_id = .data[[sample_id_column]]
        )
}

#' @rdname get_metadata
#' @keywords internal
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
#' @return A `data.frame` of detected peaks (one row per peak), or `NULL` if no
#' peaks are available.
#' @keywords internal
get_detected_peaks <- function(obj) {
    UseMethod("get_detected_peaks")
}

.get_detected_peaks_xcms <- function(obj) {
    if (is_xcms_processed_data(obj) && xcms::hasChromPeaks(obj)) {
        as.data.frame(xcms::chromPeaks(obj)) |>
            dplyr::rename(sample_index = sample)
    } else {
        NULL
    }
}

#' @rdname get_detected_peaks
#' @keywords internal
get_detected_peaks.character <- function(obj) {
    NULL
}

#' @rdname get_detected_peaks
#' @keywords internal
get_detected_peaks.XCMSnExp <- function(obj) {
    .get_detected_peaks_xcms(obj)
}

#' @rdname get_detected_peaks
#' @keywords internal
get_detected_peaks.MsExperiment <- function(obj) {
    .get_detected_peaks_xcms(obj)
}

#' @rdname get_detected_peaks
#' @keywords internal
get_detected_peaks.MChromatograms <- function(obj) {
    NULL
}

#' @rdname get_detected_peaks
#' @keywords internal
get_detected_peaks.XChromatograms <- function(obj) {
    if (any(xcms::hasChromPeaks(obj))) {
        peaks <- as.data.frame(xcms::chromPeaks(obj)) |>
            dplyr::rename(sample_index = "column") |>
            select(-dplyr::all_of("mz")) # This is already present in mz_info

        # Extract mz ranges from each row of the XChromatograms object
        mz_info <- do.call(rbind, lapply(seq_len(nrow(obj)), function(i) {
            mz_range <- MSnbase::mz(obj[i, 1L])[[1L]]
            data.frame(
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
get_detected_peaks.XChromatogram <- function(obj) {
    if (xcms::hasChromPeaks(obj)) {
        mz_range <- MSnbase::mz(obj)
        as.data.frame(xcms::chromPeaks(obj)) |>
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
get_detected_peaks.XcmsRawList <- function(obj) {
    NULL
}

#' @rdname get_detected_peaks
#' @keywords internal
get_detected_peaks.ExternalDataSource <- function(obj) {
    obj@peaks
}

#' Get the grouped peaks across samples (features) from the data object
#'
#' `get_grouped_peaks()` is an internal helper that retrieves feature-level
#' grouped peaks, i.e., chromatographic peaks aligned across samples.
#'
#' @param obj A data object containing or representing samples.
#' @return A `data.frame` of grouped (feature-level) peaks,
#' or `NULL` if not available.
#' @keywords internal
get_grouped_peaks <- function(obj) {
    UseMethod("get_grouped_peaks")
}

#' @rdname get_grouped_peaks
#' @keywords internal
get_grouped_peaks.default <- function(obj) {
    return(NULL)
}

.get_grouped_peaks_xcms <- function(obj) {
    if (is_xcms_processed_data(obj) && xcms::hasFeatures(obj)) {
        as.data.frame(xcms::featureDefinitions(obj)) |>
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
get_grouped_peaks.XCMSnExp <- function(obj) {
    .get_grouped_peaks_xcms(obj)
}

#' @rdname get_grouped_peaks
#' @keywords internal
get_grouped_peaks.XcmsExperiment <- function(obj) {
    .get_grouped_peaks_xcms(obj)
}

#' @rdname get_grouped_peaks
#' @keywords internal
get_grouped_peaks.MsExperiment <- function(obj) {
    .get_grouped_peaks_xcms(obj)
}
