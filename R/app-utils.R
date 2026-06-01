#' Compute an m/z window from a target m/z and a ppm tolerance
#'
#' Internal helper used by the Shiny app to translate a user-supplied
#' centre mass and ppm tolerance into the `mzmin`/`mzmax` pair expected by
#' [lp_chromatogram()].
#'
#' @param mz A `numeric` scalar — the target m/z.
#' @param ppm A `numeric` scalar — the symmetric tolerance in parts-per-million.
#' @return A length-2 `numeric` vector `c(mzmin, mzmax)`.
#' @keywords internal
.mz_window <- function(mz, ppm) {
    delta <- mz * ppm / 1e6
    c(mz - delta, mz + delta)
}

#' Build a features matrix for an extracted ion chromatogram
#'
#' Produces the single-row `features` matrix expected by [lp_chromatogram()]
#' from a target `mz`, a ppm tolerance, and an optional retention-time window.
#' When `rt` is `NULL` the RT columns are filled with `NA` so the full RT
#' range is plotted.
#'
#' @param mz A `numeric` scalar — the target m/z.
#' @param ppm A `numeric` scalar — the ppm tolerance.
#' @param rt A `numeric` scalar or `NULL` — the target retention time (in
#'   the same unit as the data, typically seconds).
#' @param rt_tol A `numeric` scalar — the symmetric RT tolerance. Ignored
#'   when `rt` is `NULL`.
#' @return A 1x4 `matrix` with column names `mzmin`, `mzmax`, `rtmin`, `rtmax`.
#' @keywords internal
.eic_features <- function(mz, ppm, rt = NULL, rt_tol = NULL) {
    mz_w <- .mz_window(mz, ppm)
    if (is.null(rt) || is.na(rt)) {
        rtmin <- NA_real_
        rtmax <- NA_real_
    } else {
        rtmin <- rt - rt_tol
        rtmax <- rt + rt_tol
    }
    m <- matrix(
        c(mz_w[1], mz_w[2], rtmin, rtmax),
        nrow = 1,
        dimnames = list(NULL, c("mzmin", "mzmax", "rtmin", "rtmax")))
    m
}

#' Copy uploaded Shiny files to a session temp dir under their original names
#'
#' `shiny::fileInput()` writes uploaded files to randomly-named tempfiles
#' (`/tmp/.../0.mzML`) and reports the original names in a separate `name`
#' column. The `lcmsPlot()` constructor derives sample IDs from file basenames,
#' so we copy each upload to a session-scoped directory under its original
#' name to keep sample IDs meaningful.
#'
#' @param upload A `data.frame` as produced by `shiny::fileInput()` with at
#'   least the columns `datapath` and `name`. May be `NULL`.
#' @param dir A `character` directory path to copy files into. The directory
#'   is created if it does not exist.
#' @return A `character` vector of full paths to the renamed files, in the
#'   same order as `upload`. Returns `character(0)` if `upload` is `NULL`
#'   or empty.
#' @keywords internal
.uploaded_files_to_paths <- function(upload, dir) {
    if (is.null(upload) || nrow(upload) == 0) return(character(0))
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
    out <- file.path(dir, upload$name)
    file.copy(upload$datapath, out, overwrite = TRUE)
    out
}

#' Lowercased file extension without the dot
#'
#' @param path A `character` vector of file paths.
#' @return A `character` vector of the same length, lowercased extensions.
#' @keywords internal
.file_ext <- function(path) {
    tolower(tools::file_ext(path))
}

#' Classes the Shiny app accepts when loading a `.rds` or `.RData` payload
#'
#' These match the classes `lcmsPlot()` understands as `data_obj` (see
#' `R/lcmsPlotDataContainer-class.R`). Kept as a module-level constant so
#' it's easy to extend in one place.
#'
#' @keywords internal
.SUPPORTED_DATA_CLASSES <- c(
    "XCMSnExp", "MsExperiment",
    "XChromatograms", "XChromatogram", "MChromatograms",
    "XcmsRawList", "purityA")

#' Translate a set of uploaded file paths to whatever `lcmsPlot()` expects
#'
#' Inspects file extensions to decide how each upload should be handed off
#' to [lcmsPlot()]. Pure function, no Shiny — easy to unit-test.
#'
#' Dispatch:
#' - mzML / CDF / raw paths → returned unchanged so `lcmsPlot()` can dispatch
#'   them through the existing character-vector path.
#' - A single `.cdResult` path → returned unchanged; `lcmsPlot()` opens the
#'   SQLite connection internally via [is_cd_results_path()].
#' - A single `.rds` → [readRDS()]'d and the resulting object returned.
#' - A single `.RData` / `.Rdata` / `.rda` → [load()]'d into a fresh
#'   environment; the first object inheriting from one of
#'   `.SUPPORTED_DATA_CLASSES` is returned. Errors with a clear message if
#'   no such object exists.
#'
#' Empty input, unrecognised extensions, or mixed types raise a `stop()` so
#' the caller can surface the message via [.toast_error()].
#'
#' @param paths A `character` vector of full file paths.
#' @return Either a `character` vector of paths, a `character` scalar path,
#'   or an R object — whichever is appropriate for `lcmsPlot()`.
#' @keywords internal
.load_dataset <- function(paths) {
    if (length(paths) == 0) {
        stop("No files supplied.", call. = FALSE)
    }

    exts <- .file_ext(paths)
    raw_exts <- c("mzml", "cdf", "raw")

    if (all(exts %in% raw_exts)) return(paths)

    if (length(paths) > 1) {
        stop("Only one file can be uploaded for this format. ",
             "Got: ", paste(basename(paths), collapse = ", "),
             call. = FALSE)
    }

    path <- paths[1]
    ext  <- exts[1]

    if (ext == "cdresult") return(path)

    if (ext %in% c("rds")) return(readRDS(path))

    if (ext %in% c("rdata", "rda")) {
        e <- new.env()
        load(path, envir = e)
        for (nm in ls(e, all.names = TRUE)) {
            obj <- get(nm, envir = e)
            if (any(vapply(.SUPPORTED_DATA_CLASSES,
                           function(cls) methods::is(obj, cls),
                           logical(1)))) {
                return(obj)
            }
        }
        stop("RData file does not contain any of: ",
             paste(.SUPPORTED_DATA_CLASSES, collapse = ", "), ".",
             call. = FALSE)
    }

    stop("Unsupported file type: .", ext,
         ". Supported extensions: mzML, CDF, raw, cdResult, rds, RData.",
         call. = FALSE)
}

#' Show an error toast and return invisibly
#'
#' Thin wrapper around [shinytoastr::toastr_error()] used by module servers
#' to surface `lp_*()` errors without crashing the session.
#'
#' @param msg A `character` message to display.
#' @param title A `character` title for the toast.
#' @return `NULL` invisibly.
#' @keywords internal
.toast_error <- function(msg, title = "Error") {
    if (requireNamespace("shinytoastr", quietly = TRUE)) {
        shinytoastr::toastr_error(message = msg, title = title)
    }
    invisible(NULL)
}
