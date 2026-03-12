#' A list of `xcmsRaw` objects representing multiple samples
#'
#' `XcmsRawList` is a container for one or more `xcmsRaw` objects, where each
#' element corresponds to a single sample. It is the recommended way to pass
#' in-memory raw LC-MS data to `lcmsPlot`.
#'
#' @slot data A `list` of `xcmsRaw` objects, one per sample.
#' @export
setClass(
    "XcmsRawList",
    slots = list(data = "list")
)

setValidity("XcmsRawList", function(object) {
    if (length(object@data) == 0L) {
        return("XcmsRawList must contain at least one xcmsRaw object.")
    }
    if (!all(vapply(object@data, is, logical(1L), "xcmsRaw"))) {
        return("All elements of XcmsRawList must be xcmsRaw objects.")
    }
    TRUE
})

#' Create an `XcmsRawList` object
#'
#' Wraps one or more `xcmsRaw` objects into an `XcmsRawList` container for use
#' as the `dataset` argument of `lcmsPlot`.
#'
#' @param ... One or more `xcmsRaw` objects, or a single `list` of `xcmsRaw`
#' objects.
#' @return An instance of `XcmsRawList`.
#' @export
#' @examples
#' \dontrun{
#' raw1 <- xcms::xcmsRaw("sample1.mzML")
#' raw2 <- xcms::xcmsRaw("sample2.mzML")
#' xl <- XcmsRawList(raw1, raw2)
#' }
XcmsRawList <- function(...) {
    objs <- list(...)
    # Allow passing a pre-built list: XcmsRawList(list(raw1, raw2))
    if (length(objs) == 1L &&
        is.list(objs[[1L]]) &&
        !is(objs[[1L]], "xcmsRaw")) {
        objs <- objs[[1L]]
    }
    new("XcmsRawList", data = objs)
}

#' Create an `XcmsRawList` from raw LC-MS files
#'
#' Reads one or more raw MS files using `xcms::xcmsRaw()` and wraps the
#' results in an `XcmsRawList`. Optionally runs the reads in parallel via
#' `BiocParallel`.
#'
#' @param paths A `character` vector of file paths to raw MS files
#' (e.g. `.mzML`, `.mzXML`, `.CDF`).
#' @param profstep A `numeric` value passed to `xcms::xcmsRaw()` controlling
#' the mass bin size used to build the profile matrix. Defaults to `0`, which
#' skips profile matrix generation. Only set this to a positive value if you
#' need the profile matrix for peak detection.
#' @param mslevel A `numeric` value passed to `xcms::xcmsRaw()` indicating
#' which MS level to load. `NULL` loads all levels. Defaults to `NULL`.
#' @param scanrange A length-2 `integer` vector passed to `xcms::xcmsRaw()`
#' restricting the range of scans to read. `NULL` reads all scans.
#' Defaults to `NULL`.
#' @param BPPARAM A `BiocParallelParam` object controlling parallel execution.
#' When `NULL` (default) files are read sequentially.
#' @return An `XcmsRawList` object with one `xcmsRaw` element per file.
#' @export
#' @examples
#' paths <- dir(
#'     system.file("cdf", package = "faahKO"),
#'     full.names = TRUE,
#'     recursive = TRUE
#' )[1:3]
#' xl <- create_xcms_raw_list(paths)
create_xcms_raw_list <- function(
    paths,
    profstep = 1,
    mslevel = NULL,
    scanrange = NULL,
    BPPARAM = NULL
) {
    read_one <- function(path) {
        xcms::xcmsRaw(
            path,
            profstep = profstep,
            mslevel = mslevel,
            scanrange = scanrange
        )
    }

    if (!is.null(BPPARAM)) {
        objs <- BiocParallel::bplapply(paths, read_one, BPPARAM = BPPARAM)
    } else {
        objs <- lapply(paths, read_one)
    }

    new("XcmsRawList", data = objs)
}

#' Show a summary of an `XcmsRawList` object
#'
#' @param object An instance of `XcmsRawList`.
#' @return Invisible `NULL`.
#' @export
setMethod("show", "XcmsRawList", function(object) {
    cat("Object of class XcmsRawList\n")
    cat(" Number of samples:", length(object@data), "\n")
    paths <- vapply(object@data, function(x) {
        if (length(x@filepath) > 0L) x@filepath[[1L]] else "<no path>"
    }, character(1L))
    cat(" File paths:\n")
    for (p in paths) cat("  -", p, "\n")
    invisible(NULL)
})
