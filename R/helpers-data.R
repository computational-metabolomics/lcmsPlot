#' Check whether an object is from XCMS
#'
#' @param obj The input object to check.
#' @return A `logical` value indicating whether the object is
#' an XCMS one, i.e. `XCMSnExp` or `MsExperiment`.
#' @keywords internal
is_xcms_data <- function(obj) {
    inherits(obj, c("XCMSnExp", "MsExperiment"))
}

#' Check whether an object is an msPurity purityA result
#'
#' @param obj The input object to check.
#' @return A `logical` value.
#' @keywords internal
is_mspurity_data <- function(obj) {
    inherits(obj, "purityA")
}

#' Check whether an object is an XCMS data container
#'
#' @param obj The input object to check.
#' @return A `logical` value indicating whether `obj` is an
#' XCMS experiment object.
#' @keywords internal
is_xcms_processed_data <- function(obj) {
    inherits(obj, c("XCMSnExp", "XcmsExperiment"))
}

#' Check whether a path corresponds to a Compound Discoverer results file
#'
#' @param path A `character` value giving a file system path.
#' @return A `logical` value indicating whether the path corresponds to a
#' Compound Discoverer results directory.
#' @keywords internal
is_cd_results_path <- function(path) {
    is.character(path) && length(path) == 1 && endsWith(path, ".cdResult")
}

#' Check whether an object is a Compound Discoverer database connection
#'
#' The object must inherit from `DBIConnection`, and its `dbname`
#' slot must reference a path ending in `.cdResult`.
#'
#' @param obj An object to test.
#' @return A `logical` value indicating whether the object is a Compound
#' Discoverer database connection.
#' @keywords internal
is_cd_result <- function(obj) {
    inherits(obj, "DBIConnection") && is_cd_results_path(obj@dbname)
}

#' Check whether an object is a Compound Discoverer scripting-node data source
#'
#' @param obj An object to test.
#' @return A `logical` value indicating whether the object is a
#' `CompoundDiscovererNodeSource`.
#' @keywords internal
is_cd_node_source <- function(obj) {
    inherits(obj, "CompoundDiscovererNodeSource")
}

#' Check whether an object is a LipidSearch data source
#'
#' @param obj An object to test.
#' @return A `logical` value indicating whether the object is a
#' `LipidSearchSource`.
#' @keywords internal
is_lipid_search_source <- function(obj) {
    inherits(obj, "LipidSearchSource")
}

#' Resolve the first matching column name in a data frame
#'
#' Vendor exports name the same logical field differently between versions, so
#' each field is resolved against a list of candidate column names.
#'
#' @param df A `data.frame`.
#' @param candidates A `character` vector of candidate column names.
#' @return The first candidate present in `df`, or `NA_character_` if none.
#' @keywords internal
first_matching_column <- function(df, candidates) {
    hit <- intersect(candidates, colnames(df))
    if (length(hit) > 0) hit[[1]] else NA_character_
}

#' Add a dense compound rank (by descending total area) for top-N selection
#'
#' Shared by the compound-centric data sources (Compound Discoverer scripting
#' node, LipidSearch). Compounds are ranked by the total area summed across
#' samples, so `rank == 1` is the most abundant compound.
#'
#' @param rows A `data.frame` with `name` and `into` columns.
#' @param rank_column A `character` value giving the name of the rank column to
#' add.
#' @return `rows` with the rank column appended.
#' @keywords internal
add_compound_rank <- function(rows, rank_column = "compound_rank") {
    ranks <- rows |>
        dplyr::group_by(.data$name) |>
        dplyr::summarise(
            total_into = sum(.data$into, na.rm = TRUE), .groups = "drop"
        ) |>
        dplyr::mutate(
            rank = dplyr::dense_rank(dplyr::desc(.data$total_into))
        ) |>
        dplyr::select(.data$name, .data$rank)

    names(ranks)[names(ranks) == "rank"] <- rank_column
    rows |> left_join(ranks, by = "name")
}

#' Get an `XCMSnExp` example object from the `faahKO` dataset
#'
#' @param indices A `numeric` vector of sample indices to select
#' from the `faahKO` CDF files. Defaults to the first three samples.
#' @param should_group_peaks A `logical` value indicating whether
#' to group the detected peaks.
#' @return An `XCMSnExp` object containing raw data, detected chromatographic
#' peaks, and, if requested, grouped features.
#' @export
#' @examples
#' get_XCMSnExp_object_example(indices = 1:5)
get_XCMSnExp_object_example <- function(
    indices = c(1, 2, 3),
    should_group_peaks = FALSE
) {
    cdfs <- dir(system.file("cdf", package = "faahKO"), full.names = TRUE,
                recursive = TRUE)[indices]
    sample_names <- sub(
        basename(cdfs),
        pattern = ".CDF",
        replacement = "",
        fixed = TRUE
    )

    pd <- tibble(
        sample_name = sample_names,
        sample_group = toupper(sub("[0-9]+", "", sample_names))
    )

    raw_data <- MSnbase::readMSData(
        files = cdfs,
        pdata = new("AnnotatedDataFrame", pd),
        mode = "onDisk",
        msLevel = 1
    )

    cwp <- xcms::CentWaveParam(
        peakwidth = c(20, 80),
        noise = 10000,
        prefilter = c(6, 10000)
    )
    xdata <- xcms::findChromPeaks(raw_data, param = cwp)

    if (should_group_peaks) {
        xdata <- xcms::adjustRtime(
            xdata,
            param = xcms::ObiwarpParam(binSize = 0.6))
        pdp <- xcms::PeakDensityParam(
            sampleGroups = pd$sample_group,
            minFraction = 1,
            bw = 30)
        xdata <- xcms::groupChromPeaks(xdata, param = pdp)
    }

    xdata
}

#' Merge two data frames using a row-index match
#'
#' @param a The left data frame.
#' @param b The right data frame whose row order defines the index.
#' @param index_col The name of the column in `a` that contains row indices
#' referring to `b`.
#' @return A `tibble` resulting from a left join of `a` and the indexed `b`.
#' @keywords internal
merge_by_index <- function(a, b, index_col) {
    b_mod <- b |> mutate(row_id = row_number())
    a |> left_join(b_mod, by = setNames("row_id", index_col))
}

#' Remove `NULL` elements from a list
#'
#' @param lst A `list` from which `NULL` entries should be removed.
#' @return A `list` containing only the non-`NULL` elements of `lst`.
#' @keywords internal
remove_null_elements <- function(lst) {
    return(Filter(Negate(is.null), lst))
}

#' Compute an m/z range given a ppm tolerance
#'
#' @param mz A `numeric` value indicating the m/z value
#' to calculate the range for.
#' @param ppm A `numeric` value indicating the tolerance
#' in parts per million (ppm). Default is 5.
#' @return A `numeric` vector with two values indicating the m/z range.
#' @keywords internal
get_mz_range <- function(mz, ppm = 5) {
    mzdev <- mz * (ppm / 1000000)
    return(c(mz - mzdev, mz + mzdev))
}

#' Get a feature's m/z and RT ranges given different feature specifications
#'
#' @param feature The input feature which can be a vector or a data frame row.
#' The accepted columns combinations are:
#' - mz, rt (optional)
#' - mzmin, mzmax, rtmin (optional), rtmax (optional)
#' @param options The plot object's options.
#' @param full_rt_range The full RT range if an RT range is not given.
#' @return A named list defining a feature with names: feature_id, mzr, rtr.
#' @keywords internal
get_feature_data <- function(feature, options, full_rt_range) {
    # Helper: safely extract a value by name from vector or tibble
    get_val <- function(x, name) {
        if (is.data.frame(x)) {
            if (name %in% names(x)) return(x[[name]][1])
        } else if (is.vector(x)) {
            if (name %in% names(x)) return(x[[name]])
        }
        return(NA)
    }

    mz <- get_val(feature, "mz")
    rt <- get_val(feature, "rt")
    mzmin  <- get_val(feature, "mzmin")
    mzmax  <- get_val(feature, "mzmax")
    rtmin  <- get_val(feature, "rtmin")
    rtmax  <- get_val(feature, "rtmax")

    # Case 1: Single mz (± tolerance) with optional rt
    if (!is.na(mz)) {
        mzr <- get_mz_range(mz, options$chromatograms$ppm)

        if (!is.na(rt)) {
            rtr <- c(
                rt - options$chromatograms$rt_tol,
                rt + options$chromatograms$rt_tol)
            feature_id <- sprintf("M%dT%d", round(mz), round(rt))
        } else {
            rtr <- full_rt_range
            feature_id <- sprintf("M%d", round(mz))
        }

        # Case 2: Explicit mz range, with optional rt range
    } else if (!is.na(mzmin) && !is.na(mzmax)) {
        mzr <- c(mzmin, mzmax)

        if (!is.na(rtmin) && !is.na(rtmax)) {
            rtr <- c(rtmin, rtmax)
            feature_id <- sprintf(
                "M%dT%d",
                round(mean(c(mzmin, mzmax))),
                round(mean(c(rtmin, rtmax))))
        } else {
            rtr <- full_rt_range
            feature_id <- sprintf("M%d", round(mean(c(mzmin, mzmax))))
        }

    } else {
        stop("Unsupported feature format: must provide either 'mz' or ('mzmin' and 'mzmax').")
    }

    list(
        feature_id = feature_id,
        mzr = mzr,
        rtr = rtr
    )
}

#' Get the feature data (m/z and RT ranges) for a set of features
#'
#' @param options The plot object's options. This contains the input features.
#' @param sample_metadata The sample's metadata.
#' @param grouped_peaks The grouped peaks to use as input features.
#' @param full_rt_range The full RT range if an RT range is not given.
#' @return A list of feature data (see get_feature_data).
#' @keywords internal
get_features <- function(
    options,
    sample_metadata,
    grouped_peaks = NULL,
    full_rt_range = NULL
) {
    if (is.null(grouped_peaks)) {
        input_features <- options$chromatograms$features
    } else {
        input_features <- grouped_peaks
    }

    if (is.data.frame(input_features) &&
        "sample_id" %in% colnames(input_features)) {
        feature_indices <- which(
            input_features$sample_id == sample_metadata$sample_id)
    } else {
        feature_indices <- seq_len(nrow(input_features))
    }

    features <- list()

    for (i in seq_len(length(feature_indices))) {
        feature_index <- feature_indices[i]
        feature <- input_features[feature_index, ]
        features[[i]] <- get_feature_data(feature, options, full_rt_range)
    }

    return(features)
}

#' Retrieve raw and adjusted retention times from an xcms object
#'
#' Extracts raw and adjusted retention times from an xcms-processed object
#' when retention time correction has been performed.
#'
#' If the object is not an xcms processed data object, or if retention time
#' adjustment has not been applied, the function returns \code{NULL}.
#'
#' @param obj An object potentially containing xcms-processed LC-MS data.
#' @return A `tibble` with one row per scan, containing:
#' \describe{
#'   \item{file_index}{Index of the originating raw data file.}
#'   \item{raw_rt}{Original (unadjusted) retention time.}
#'   \item{adj_rt}{Adjusted retention time after RT correction.}
#' }
#' If no adjusted retention times are available, \code{NULL} is returned.
#' @keywords internal
get_adjusted_rts <- function(obj) {
    if (is_xcms_processed_data(obj) && xcms::hasAdjustedRtime(obj)) {
        adjusted_rt <- xcms::adjustedRtime(obj)
        file_indices <- xcms::fromFile(obj)
        return(tibble(
            file_index = file_indices,
            raw_rt = xcms::rtime(obj, adjusted = FALSE),
            adj_rt = adjusted_rt
        ))
    } else {
        return(NULL)
    }
}

#' Detect field separator from file extension
#'
#' Determines the column separator to use when reading a delimited text file
#' based on its file extension. Files ending in `.tsv` (case-insensitive) are
#' assumed to be tab-delimited; all others default to comma-delimited.
#'
#' @param path A `character` value indicating the path to the file whose
#' separator should be detected.
#' @return A `character` value indicating the field separator: `"\t"` for TSV
#' files, otherwise `","`.
#' @keywords internal
detect_separator <- function(path) {
    if (grepl("\\.tsv$", path, ignore.case = TRUE)) "\t" else ","
}
