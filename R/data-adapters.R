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

    xcms::phenoData(obj)@data |>
        mutate(
            sample_index = row_number(),
            sample_id = .data[[sample_id_column]],
            sample_path = xcms::fileNames(obj)
        )
}

#' @rdname get_metadata
#' @keywords internal
get_metadata.MsExperiment <- function(obj, sample_id_column, metadata) {
    if (!is.null(metadata)) {
        MsExperiment::sampleData(obj) <- metadata
    }

    MsExperiment::sampleData(obj) |>
        as.data.frame() |>
        mutate(
            sample_index = row_number(),
            sample_id = .data[[sample_id_column]],
            sample_path = xcms::fileNames(obj)
        )
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
            sample_id = .data$StudyFileID
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
    return(NULL)
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
