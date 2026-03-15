#' Virtual base class for raw MS file readers
#'
#' `MsRawReader` is an abstract S4 class that defines the interface for
#' reading raw mass spectrometry data files. Concrete subclasses implement
#' the interface for specific backends (e.g., `mzR`, `rawrr`).
#'
#' @keywords internal
setClass("MsRawReader", contains = "VIRTUAL")

#' An `mzR`-backed raw MS file reader
#'
#' Wraps an open `mzR` connection (as returned by `mzR::openMSfile()`).
#'
#' @slot connection An `mzR` connection object.
#' @keywords internal
setClass(
    "MzrReader",
    contains = "MsRawReader",
    slots = list(connection = "ANY")
)

#' A `rawrr`-backed raw MS file reader
#'
#' Stores the path to a ThermoFisher `.raw` file. `rawrr` does not maintain
#' persistent connections, so the path is re-used for each read operation.
#'
#' @slot path A `character` value giving the path to the `.raw` file.
#' @keywords internal
setClass(
    "RawrrReader",
    contains = "MsRawReader",
    slots = list(path = "character")
)

#' Get scan header information from a raw MS file reader
#'
#' @param reader An instance of `MsRawReader`.
#' @return A `tibble` with at least the columns `seqNum`,
#' `retentionTime`, `msLevel`, `basePeakIntensity`, and `totIonCurrent`.
#' @keywords internal
setGeneric("ms_header", function(reader) standardGeneric("ms_header"))

#' Get peak matrices for selected scans from a raw MS file reader
#'
#' @param reader An instance of `MsRawReader`.
#' @param scans An `integer` vector of scan sequence numbers to retrieve.
#' @return A `list` of two-column matrices with columns `mz` and `intensity`,
#' one element per scan (always a list, even for a single scan).
#' @keywords internal
setGeneric(
    "ms_peaks",
    function(reader, scans) standardGeneric("ms_peaks")
)

#' Close a raw MS file reader connection
#'
#' @param reader An instance of `MsRawReader`.
#' @return `invisible(NULL)`
#' @keywords internal
setGeneric("ms_close", function(reader) standardGeneric("ms_close"))

#' Extract an extracted ion chromatogram (XIC) using the native backend
#'
#' `ms_chromatogram()` uses the backend's own XIC extraction mechanism rather
#' than reconstructing chromatograms scan-by-scan from `ms_peaks()`.
#' Currently only implemented for `RawrrReader` (ThermoFisher `.raw` files),
#' which delegates to `rawrr::readChromatogram()`.
#'
#' @param reader An instance of `MsRawReader`.
#' @param mz A `numeric` scalar — the target m/z.
#' @param ppm A `numeric` scalar — the mass tolerance in ppm.
#' @param rt_range A length-2 `numeric` vector — the RT window in seconds.
#' @return A `list` with two `tibble`s: `chromatograms` (columns `rt` and
#' `intensity`) and `mass_traces` (empty for backends that do not expose
#' per-scan m/z data via this interface).
#' @keywords internal
setGeneric(
    "ms_chromatogram",
    function(reader, mz, ppm, rt_range) standardGeneric("ms_chromatogram")
)

#' @rdname ms_header
setMethod("ms_header", "MzrReader", function(reader) {
    hdr <- mzR::header(reader@connection)

    if (is(reader@connection, "mzRnetCDF")) {
        ms1_idx <- which(hdr$msLevel == 1)
        scan_nums <- hdr$seqNum[ms1_idx]
        num_scans <- length(scan_nums)

        bpi <- numeric(num_scans)
        tic_vals <- numeric(num_scans)

        batch_size <- 100
        batch_starts <- seq(1, num_scans, by = batch_size)

        for (start in batch_starts) {
            end <- min(start + batch_size - 1, num_scans)
            batch_scans <- scan_nums[start:end]
            peaks_list <- mzR::peaks(reader@connection, scans = batch_scans)

            for (i in seq_along(batch_scans)) {
                pk <- peaks_list[[i]]
                has_data <- length(pk) > 0 && nrow(pk) > 0
                bpi[start + i - 1] <- if (has_data) max(pk[, 2]) else 0
                tic_vals[start + i - 1] <- if (has_data) sum(pk[, 2]) else 0
            }
        }

        hdr$basePeakIntensity[ms1_idx] <- bpi
        hdr$totIonCurrent[ms1_idx] <- tic_vals
    }

    hdr
})

#' @rdname ms_peaks
setMethod("ms_peaks", "MzrReader", function(reader, scans) {
    result <- mzR::peaks(reader@connection, scans = scans)
    if (is.matrix(result)) list(result) else result
})

#' @rdname ms_chromatogram
setMethod(
    "ms_chromatogram", "MzrReader",
    function(reader, mz, ppm, rt_range) {
        stop("ms_chromatogram() is not implemented for MzrReader.")
    }
)

#' @rdname ms_close
setMethod("ms_close", "MzrReader", function(reader) {
    mzR::close(reader@connection)
    invisible(NULL)
})

#' @rdname ms_header
setMethod("ms_header", "RawrrReader", function(reader) {
    .rawrr_ms_order_to_level <- function(ms_order) {
        ifelse(ms_order == "Ms", 1L, as.integer(sub("Ms", "", ms_order)))
    }

    idx <- rawrr::readIndex(reader@path)
    bpc <- rawrr::readChromatogram(rawfile = reader@path, type = "bpc")
    tic <- rawrr::readChromatogram(rawfile = reader@path, type = "tic")
    tibble(
        seqNum = idx$scan,
        retentionTime = idx$StartTime * 60,
        msLevel = .rawrr_ms_order_to_level(idx$MSOrder),
        basePeakIntensity = bpc$intensities,
        totIonCurrent = tic$intensities
    )
})

#' @rdname ms_peaks
setMethod("ms_peaks", "RawrrReader", function(reader, scans) {
    spectra <- rawrr::readSpectrum(reader@path, scan = scans)
    lapply(spectra, function(s) cbind(mz = s$mZ, intensity = s$intensity))
})

#' @rdname ms_chromatogram
setMethod(
    "ms_chromatogram", "RawrrReader",
    function(reader, mz, ppm, rt_range) {
        chrom <- rawrr::readChromatogram(
            reader@path, type = "xic", mass = mz, tol = ppm)
        rt_seconds <- chrom[[1]]$times * 60
        in_range <- rt_seconds >= rt_range[1] & rt_seconds <= rt_range[2]
        tibble(
            rt = rt_seconds[in_range],
            intensity = chrom[[1]]$intensities[in_range]
        )
    }
)

#' @rdname ms_close
setMethod("ms_close", "RawrrReader", function(reader) {
    invisible(NULL)
})

#' An `xcmsRaw`-backed raw MS file reader
#'
#' Wraps an in-memory `xcmsRaw` object. No file connection is opened or closed.
#'
#' @slot obj An `xcmsRaw` object.
#' @keywords internal
setClass(
    "XcmsRawReader",
    contains = "MsRawReader",
    slots = list(obj = "ANY")
)

#' @rdname ms_header
setMethod("ms_header", "XcmsRawReader", function(reader) {
    obj <- reader@obj
    scanidx <- obj@scanindex
    ints <- obj@env$intensity
    nscans <- length(scanidx)

    bpi <- numeric(nscans)
    for (i in seq_len(nscans)) {
        start <- scanidx[i] + 1L
        end <- if (i < nscans) scanidx[i + 1L] else length(ints)
        bpi[i] <- if (start <= end && length(ints) > 0) max(ints[start:end]) else 0
    }

    tibble(
        seqNum = seq_len(nscans),
        retentionTime = obj@scantime,
        msLevel = 1L,
        basePeakIntensity = bpi,
        totIonCurrent = obj@tic
    )
})

#' @rdname ms_peaks
setMethod("ms_peaks", "XcmsRawReader", function(reader, scans) {
    obj <- reader@obj
    scanidx <- obj@scanindex
    mzs <- obj@env$mz
    ints <- obj@env$intensity
    nscans <- length(scanidx)

    lapply(scans, function(s) {
        start <- scanidx[s] + 1L
        end <- if (s < nscans) scanidx[s + 1L] else length(mzs)
        if (start <= end) {
            cbind(mz = mzs[start:end], intensity = ints[start:end])
        } else {
            matrix(
                numeric(0),
                ncol = 2,
                dimnames = list(NULL, c("mz", "intensity"))
            )
        }
    })
})

#' @rdname ms_chromatogram
setMethod(
    "ms_chromatogram", "XcmsRawReader",
    function(reader, mz, ppm, rt_range) {
        stop("ms_chromatogram() is not implemented for XcmsRawReader.")
    }
)

#' @rdname ms_close
setMethod("ms_close", "XcmsRawReader", function(reader) {
    invisible(NULL)
})

#' Convert an `XcmsRawList` to a named list of `XcmsRawReader`s
#'
#' Each element is keyed by the file path stored in `@filepath` of the
#' corresponding `xcmsRaw` object.
#'
#' @param xcmsraw_list An `XcmsRawList` object.
#' @return A named `list` of `XcmsRawReader` objects.
#' @keywords internal
xcmsraw_to_readers <- function(xcmsraw_list) {
    objs <- xcmsraw_list@data
    paths <- vapply(objs, function(obj) {
        if (length(obj@filepath) > 0L) obj@filepath[[1L]] else ""
    }, character(1L))
    readers <- lapply(objs, function(obj) new("XcmsRawReader", obj = obj))
    setNames(readers, paths)
}

#' Open a raw MS file as an `MsRawReader`
#'
#' Dispatches to `MzrReader` for standard open formats (`.mzML`, `.mzXML`,
#' `.CDF`) or `RawrrReader` for ThermoFisher `.raw` files.
#'
#' @param path A `character` value giving the file path.
#' @return An instance of `MsRawReader` (either `MzrReader` or `RawrrReader`).
#' @keywords internal
open_raw_reader <- function(path) {
    ext <- tolower(tools::file_ext(path))
    if (ext == "raw") {
        if (!requireNamespace("rawrr", quietly = TRUE)) {
            stop(
                "Package 'rawrr' is required to read ThermoFisher .raw files. ",
                "Install it with: BiocManager::install('rawrr')"
            )
        }
        new("RawrrReader", path = path)
    } else {
        new("MzrReader", connection = mzR::openMSfile(path))
    }
}

#' Get `MsRawReader` objects for a set of sample paths
#'
#' @param sample_paths A `character` vector of file paths to the
#' raw MS data files (e.g., `.mzML`, `.mzXML`, `.CDF`, `.raw`).
#' @return A named list of `MsRawReader` objects, where each element
#' corresponds to a file in `sample_paths`.
#' @keywords internal
io_get_raw_data <- function(sample_paths) {
    raw_data <- list()

    for (sample_path in sample_paths) {
        raw_data[[sample_path]] <- open_raw_reader(sample_path)
    }

    return(raw_data)
}

#' Close open `MsRawReader` connections
#'
#' @param raw_data A list of `MsRawReader` objects (as returned by
#' `io_get_raw_data()`).
#' @return `NULL`
#' @keywords internal
io_close_raw_data <- function(raw_data) {
    for (raw_obj in raw_data) {
        ms_close(raw_obj)
    }
}
