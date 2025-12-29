compounds_results_columns_map <- c(
    sample_id = "StudyFileID",
    id = "ID",
    name = "Name",
    formula = "ElementalCompositionFormula",
    adduct = "ReferenceIon",
    trace = "Trace",
    rt = "Ion_RT",
    rtmin = "LeftRT",
    rtmax = "RightRT",
    mz = "MassOverCharge",
    maxo = "Intensity",
    into = "Area"
)

#' Open a connection to a Compound Discoverer results database
#'
#' Establishes a read-only SQLite database connection to a Compound Discoverer
#' results directory.
#'
#' @param cd_result_path A `character` value giving the path to a
#' Compound Discoverer results file.
#' @return A `DBIConnection` object connected to the Compound Discoverer
#' SQLite database.
#' @keywords internal
open_cd_result_connection <- function(cd_result_path) {
    DBI::dbConnect(RSQLite::SQLite(), cd_result_path)
}

#' Retrieve workflow input files from a Compound Discoverer database
#'
#' Extracts the contents of the `WorkflowInputFiles` table from a
#' Compound Discoverer results database.
#'
#' @param conn A `DBIConnection` to a Compound Discoverer results database.
#' @return A `data.frame` containing information about the workflow input files.
#' @keywords internal
get_workflow_input_files <- function(conn) {
    tbl(conn, "WorkflowInputFiles") |>
        collect()
}

#' Extract XIC traces associated with selected compounds
#'
#' Queries a Compound Discoverer results database to retrieve extracted ion
#' chromatogram (XIC) traces associated with consolidated compounds matching
#' a user-supplied filter expression.
#'
#' The query joins multiple internal Compound Discoverer tables to link
#' consolidated compounds to reference ions, chromatogram peaks, and XIC
#' trace data.
#'
#' @param conn A `DBIConnection` to a Compound Discoverer results database.
#' @param compounds_query_str A `character` value giving a filtering expression
#' evaluated on the resulting compound table (e.g. using compound name,
#' formula, retention time, or m/z).
#' @return A `data.frame` containing XIC trace metadata and binary trace data
#' for the selected compounds.
#' @keywords internal
get_xic_traces_from_compounds <- function(conn, compounds_query_str) {
    compounds_query <- rlang::parse_expr(compounds_query_str)

    # 1. Define the table references
    ion_tab <- tbl(conn, "UnknownCompoundIonInstanceItems") |>
        select(
            IonID = .data$ID,
            .data$IonDescription,
            Ion_RT = .data$RetentionTime,
            .data$Intensity,
            .data$Area,
            .data$StudyFileID)

    inst_tab <- tbl(conn, "UnknownCompoundInstanceItems") |>
        select(InstanceID = .data$ID, .data$ReferenceIon)

    comp_tab <- tbl(conn, "ConsolidatedUnknownCompoundItems") |>
        select(
            CompoundID = .data$ID,
            .data$Name,
            .data$ElementalCompositionFormula,
            .data$MassOverCharge)

    peak_tab  <- tbl(conn, "ChromatogramPeakItems") |>
        select(
            PeakID = .data$ID,
            .data$LeftRT,
            .data$RightRT,
            .data$IsRefPeak)

    # 2. Define link tables
    ion_xic_link <- tbl(
        conn, "UnknownCompoundIonInstanceItemsXicTraceItems")
    inst_ion_link <- tbl(
        conn, "UnknownCompoundInstanceItemsUnknownCompoundIonInstanceItems")
    cons_inst_link <- tbl(
        conn, "ConsolidatedUnknownCompoundItemsUnknownCompoundInstanceItems")
    ion_peak_link <- tbl(
        conn, "UnknownCompoundIonInstanceItemsChromatogramPeakItems")

    xic_items <- tbl(conn, "XicTraceItems")

    # 3. Build the pipeline
    query_result <- xic_items |>
        # Join XIC to Ion via Link
        inner_join(
            ion_xic_link,
            by = c("ID" = "XicTraceItemsID")) |>
        inner_join(
            ion_tab,
            by = c("UnknownCompoundIonInstanceItemsID" = "IonID")) |>

        # Join Ion to Chromatogram Peaks
        inner_join(
            ion_peak_link,
            by = "UnknownCompoundIonInstanceItemsID") |>
        inner_join(
            peak_tab,
            by = c("ChromatogramPeakItemsID" = "PeakID")) |>

        # Join Ion to Instance via Link
        inner_join(
            inst_ion_link,
            by = "UnknownCompoundIonInstanceItemsID") |>
        inner_join(
            inst_tab,
            by = c("UnknownCompoundInstanceItemsID" = "InstanceID")) |>

        # Join Instance to Consolidated via Link
        inner_join(
            cons_inst_link,
            by = "UnknownCompoundInstanceItemsID") |>
        inner_join(
            comp_tab,
            by = c("ConsolidatedUnknownCompoundItemsID" = "CompoundID")) |>

        # Filters
        filter(.data$IonDescription == .data$ReferenceIon) |>
        filter(.data$IsRefPeak == TRUE) |>

        # Select and rename columns
        select(
            .data$StudyFileID,
            .data$ID,
            .data$Name,
            .data$ElementalCompositionFormula,
            .data$ReferenceIon,
            .data$Trace,
            .data$Ion_RT,
            .data$LeftRT,
            .data$RightRT,
            .data$MassOverCharge,
            .data$Intensity,
            .data$Area
        ) |>
        rename(!!!compounds_results_columns_map) |>
        mutate(
            rt = .data$rt * 60,
            rtmin = .data$rtmin * 60,
            rtmax = .data$rtmax * 60
        ) |>
        filter(!!compounds_query)

    # 4. Collect
    query_result |> collect()
}

#' Parse a binary XIC trace from Compound Discoverer
#'
#' Decodes a binary XIC trace blob stored in a Compound Discoverer results
#' database and converts it into a retention time-intensity data frame.
#'
#' @param data A raw vector or binary object containing the encoded XIC trace.
#' @return A `data.frame` with columns:
#' \describe{
#'   \item{rt}{Retention time in seconds.}
#'   \item{intensity}{Signal intensity.}
#' }
#' If the trace cannot be parsed, \code{NULL} is returned.
#' @keywords internal
parse_trace <- function(data) {
    # check data
    if (is.null(data) || length(data) == 0) {
        return(NULL)
    }

    barray <- as.raw(data)

    # check for old XML format (PK)
    if (barray[1] == as.raw(0x50) && barray[2] == as.raw(0x4b)) {
        message("Old XML trace format is not supported!")
        return(NULL)
    }

    # unzip data (gzip)
    if (barray[1] == as.raw(0x1f) && barray[2] == as.raw(0x8b)) {
        barray <- memDecompress(barray, type = "gzip")
    }

    # remove leading zero byte
    if (barray[1] == as.raw(0x00)) {
        barray <- barray[-1]
    }

    # get version (byte 17 in Python -> index 17 in R)
    version <- as.integer(barray[17])

    # get converted GUID
    guid_bytes <- c(
        rev(barray[seq_len(4)]),
        rev(barray[4 + seq_len(2)]),
        rev(barray[6 + seq_len(2)]),
        barray[8 + seq_len(9)]
    )

    guid_hex <- toupper(sprintf("%02X", as.integer(guid_bytes)))
    guid <- paste0(
        paste0(guid_hex[seq_len(8)], collapse = ""),
        "-",
        paste0(guid_hex[8 + seq_len(4)], collapse = ""),
        "-",
        paste0(guid_hex[12 + seq_len(4)], collapse = ""),
        "-",
        paste0(guid_hex[16 + seq_len(4)], collapse = ""),
        "-",
        paste0(guid_hex[20 + seq_len(12)], collapse = "")
    )

    # move to trace data
    barray <- barray[-seq_len(17)]
    if (length(barray) == 0) {
        return(NULL)
    }

    # number of points
    count <- readBin(barray[seq_len(4)], integer(), size = 4, endian = "little")

    # initialize buffers
    spectra     <- rep(NA_integer_, count)
    times       <- numeric(0)
    intensities <- numeric(0)
    noise       <- rep(NA_real_, count)
    masses      <- rep(NA_real_, count)

    i <- 5  # R index after reading count

    # spectrum IDs
    if (barray[i] > 0) {
        size <- 4 * count
        spectra <- readBin(
            barray[(i + seq_len(size))],
            integer(),
            size = 4,
            endian = "little",
            n = count
        )
        i <- i + size
    }
    i <- i + 1

    # RTs (minutes)
    if (barray[i] > 0) {
        size <- 4 * count
        times <- readBin(
            barray[(i + seq_len(size))],
            numeric(),
            size = 4,
            endian = "little",
            n = count
        )
        i <- i + size
    }
    i <- i + 1

    # intensities
    if (barray[i] > 0) {
        size <- 4 * count
        intensities <- readBin(
            barray[(i + seq_len(size))],
            numeric(),
            size = 4,
            endian = "little",
            n = count
        )
        i <- i + size
    }
    i <- i + 1

    # noise
    if (barray[i] > 0) {
        size <- 4 * count
        noise <- readBin(
            barray[(i + seq_len(size))],
            numeric(),
            size = 4,
            endian = "little",
            n = count
        )
        i <- i + size
    }
    i <- i + 1

    # masses (version > 1)
    if (version > 1 && barray[i] > 0) {
        size <- 4 * count
        masses <- readBin(
            barray[(i + seq_len(size))],
            numeric(),
            size = 4,
            endian = "little",
            n = count
        )
    }

    # build points
    data.frame(
        rt = times * 60,
        intensity = intensities
    )
}
