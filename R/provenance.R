#' The provenance of a figure or report
#'
#' Records what an output of lcmsPlot was made from, so that it can be traced
#' back and reproduced: the versions of lcmsPlot, R, Bioconductor and every
#' loaded package, the layers of the pipeline with their arguments (from the
#' `history` slot), and the input files with their size and modification time.
#'
#' @param object An instance of class `lcmsPlotClass`.
#' @param checksums A `logical` value indicating whether to record the MD5
#' checksum of every input file. Reading large raw files can take a while.
#' @param output A named `list` describing the output (file, format, size).
#' @return A named `list`, ready to be written as JSON.
#' @keywords internal
.provenance <- function(object, checksums = FALSE, output = list()) {
    options <- object@options
    list(
        generator = paste("lcmsPlot", .package_version("lcmsPlot")),
        generated_at = format(
            Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
        output = output,
        environment = list(
            lcmsPlot = .package_version("lcmsPlot"),
            r = R.version$version.string,
            platform = R.version$platform,
            bioconductor = .package_version("BiocVersion"),
            packages = .loaded_package_versions()),
        pipeline = list(
            sample_id_column = options$sample_id_column,
            batch_size = options$batch_size,
            batch_index = options$batch_index,
            layers = lapply(object@history, function(item) {
                list(layer = item$name, args = .json_safe(item$args))
            })),
        input = list(
            source = class(object@data@source)[1],
            files = .input_files(sampleMetadata(object), checksums)))
}

.package_version <- function(package) {
    tryCatch(
        as.character(utils::packageVersion(package)),
        error = function(e) NA_character_)
}

# The version of every attached and loaded package, by name.
.loaded_package_versions <- function() {
    info <- utils::sessionInfo()
    packages <- c(info$otherPkgs, info$loadedOnly)
    versions <- vapply(packages, `[[`, character(1), "Version")
    as.list(versions[order(names(versions))])
}

# The input files of the samples, with their size and modification time, and
# their MD5 checksum when `checksums` is TRUE. Files that do not exist on disk
# (such as samples of an in-memory object) are listed with NA.
.input_files <- function(samples, checksums = FALSE) {
    if (!"sample_path" %in% colnames(samples)) {
        return(list())
    }
    paths <- as.character(samples$sample_path)
    exists <- !is.na(paths) & file.exists(paths)
    info <- file.info(paths[exists])

    files <- data.frame(
        sample_id = as.character(samples$sample_id),
        path = paths,
        size = NA_real_,
        modified = NA_character_)
    files$size[exists] <- info$size
    files$modified[exists] <- format(
        info$mtime, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    if (isTRUE(checksums)) {
        files$md5 <- NA_character_
        files$md5[exists] <- unname(tools::md5sum(paths[exists]))
    }
    files
}

# A value that jsonlite can write: functions and formulas as their code,
# matrices as data frames, and other objects by their class.
.json_safe <- function(x) {
    if (is.null(x)) {
        NULL
    } else if (is.function(x) || inherits(x, "formula")) {
        paste(deparse(x), collapse = "\n")
    } else if (is.data.frame(x) || is(x, "DataFrame")) {
        df <- as.data.frame(x)
        df[] <- lapply(df, function(col) {
            if (is.factor(col)) as.character(col) else col
        })
        df
    } else if (is.matrix(x) && !is.null(colnames(x))) {
        as.data.frame(x)
    } else if (is.factor(x)) {
        as.character(x)
    } else if (inherits(x, c("Date", "POSIXt"))) {
        format(x)
    } else if (is.list(x) && !is.object(x)) {
        lapply(x, .json_safe)
    } else if (is.atomic(x) && !is.object(x)) {
        x
    } else {
        list(class = class(x)[1])
    }
}

# The provenance file written next to an output: `figure.pdf` gives
# `figure.pdf.provenance.json`, so that the PDF and PNG of a figure each keep
# their own.
.provenance_path <- function(file) {
    paste0(file, ".provenance.json")
}

.provenance_json <- function(prov, pretty = TRUE) {
    jsonlite::toJSON(
        prov, auto_unbox = TRUE, pretty = pretty, null = "null",
        na = "null", digits = NA)
}

.write_provenance <- function(prov, file) {
    path <- .provenance_path(file)
    writeLines(.provenance_json(prov), path)
    invisible(path)
}
