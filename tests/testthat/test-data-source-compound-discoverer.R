# Build a minimal Compound Discoverer .cdResult SQLite database containing only
# the tables joined by get_xic_traces_from_compounds(), and return an open
# connection to it.
#
# Two compounds, each detected in two study files. `checked` controls whether
# the Compounds table carries the `Checked` column at all: Compound Discoverer
# only materialises it once compounds have been checked and the result saved in
# the application, so both variants occur in the wild.
.write_cd_result_db <- function(checked = c(1L, 0L)) {
    path <- tempfile(fileext = ".cdResult")
    conn <- DBI::dbConnect(RSQLite::SQLite(), path)

    compounds <- data.frame(
        ID = c(1L, 2L),
        Name = c("Betaine", "Creatine"),
        ElementalCompositionFormula = c("C5 H11 N O2", "C4 H9 N3 O2"),
        MassOverCharge = c(118.0863, 132.0766)
    )
    if (!is.null(checked)) {
        compounds$Checked <- checked
    }
    DBI::dbWriteTable(conn, "ConsolidatedUnknownCompoundItems", compounds)

    # One instance per compound, both referencing the [M+H]+1 ion.
    DBI::dbWriteTable(conn, "UnknownCompoundInstanceItems", data.frame(
        ID = c(1L, 2L),
        ReferenceIon = c("[M+H]+1", "[M+H]+1")
    ))

    # Each compound is detected in both study files.
    DBI::dbWriteTable(conn, "UnknownCompoundIonInstanceItems", data.frame(
        ID = seq_len(4),
        IonDescription = rep("[M+H]+1", 4),
        RetentionTime = c(0.86, 0.86, 0.90, 0.90),
        Intensity = c(500, 300, 800, 600),
        Area = c(5000, 3000, 8000, 6000),
        StudyFileID = c(1L, 2L, 1L, 2L)
    ))

    DBI::dbWriteTable(conn, "ChromatogramPeakItems", data.frame(
        ID = seq_len(4),
        LeftRT = c(0.80, 0.80, 0.84, 0.84),
        RightRT = c(0.92, 0.92, 0.96, 0.96),
        IsRefPeak = rep(1L, 4)
    ))

    # Trace blobs are never decoded here, so a placeholder is enough.
    DBI::dbWriteTable(conn, "XicTraceItems", data.frame(
        ID = seq_len(4),
        Trace = I(replicate(4, as.raw(c(0x00, 0x01)), simplify = FALSE))
    ))

    DBI::dbWriteTable(
        conn, "UnknownCompoundIonInstanceItemsXicTraceItems", data.frame(
            UnknownCompoundIonInstanceItemsID = seq_len(4),
            XicTraceItemsID = seq_len(4)
        ))

    DBI::dbWriteTable(
        conn, "UnknownCompoundIonInstanceItemsChromatogramPeakItems",
        data.frame(
            UnknownCompoundIonInstanceItemsID = seq_len(4),
            ChromatogramPeakItemsID = seq_len(4)
        ))

    DBI::dbWriteTable(
        conn, "UnknownCompoundInstanceItemsUnknownCompoundIonInstanceItems",
        data.frame(
            UnknownCompoundInstanceItemsID = c(1L, 1L, 2L, 2L),
            UnknownCompoundIonInstanceItemsID = seq_len(4)
        ))

    DBI::dbWriteTable(
        conn, "ConsolidatedUnknownCompoundItemsUnknownCompoundInstanceItems",
        data.frame(
            ConsolidatedUnknownCompoundItemsID = c(1L, 2L),
            UnknownCompoundInstanceItemsID = c(1L, 2L)
        ))

    conn
}

test_that("has_compound_checked_column() detects the Checked column", {
    conn <- .write_cd_result_db(checked = c(1L, 0L))
    on.exit(DBI::dbDisconnect(conn), add = TRUE)

    expect_true(has_compound_checked_column(conn))

    conn_none <- .write_cd_result_db(checked = NULL)
    on.exit(DBI::dbDisconnect(conn_none), add = TRUE)

    expect_false(has_compound_checked_column(conn_none))
})

test_that("compounds_query can filter on the checked column", {
    conn <- .write_cd_result_db(checked = c(1L, 0L))
    on.exit(DBI::dbDisconnect(conn), add = TRUE)

    res <- get_xic_traces_from_compounds(conn, "checked")

    expect_type(res$checked, "logical")
    expect_true(all(res$checked))
    expect_equal(sort(unique(res$name)), "Betaine")
    # Betaine is detected in both study files.
    expect_equal(nrow(res), 2)
})

test_that("compounds_query can filter on the negated checked column", {
    conn <- .write_cd_result_db(checked = c(1L, 0L))
    on.exit(DBI::dbDisconnect(conn), add = TRUE)

    res <- get_xic_traces_from_compounds(conn, "!checked")

    expect_false(any(res$checked))
    expect_equal(sort(unique(res$name)), "Creatine")
})

test_that("checked combines with other compound columns", {
    conn <- .write_cd_result_db(checked = c(1L, 1L))
    on.exit(DBI::dbDisconnect(conn), add = TRUE)

    res <- get_xic_traces_from_compounds(conn, "checked & into > 5000")

    expect_equal(sort(unique(res$name)), "Creatine")
    expect_true(all(res$into > 5000))
})

test_that("referencing checked without the column raises a clear error", {
    conn <- .write_cd_result_db(checked = NULL)
    on.exit(DBI::dbDisconnect(conn), add = TRUE)

    expect_error(
        get_xic_traces_from_compounds(conn, "checked"),
        "'checked' is not available in this .cdResult file",
        fixed = TRUE
    )
})

test_that("queries not referencing checked work without the column", {
    conn <- .write_cd_result_db(checked = NULL)
    on.exit(DBI::dbDisconnect(conn), add = TRUE)

    res <- get_xic_traces_from_compounds(conn, 'name == "Betaine"')

    expect_equal(unique(res$name), "Betaine")
    expect_false("checked" %in% colnames(res))
})

test_that("a NULL compounds_query applies no filter", {
    conn <- .write_cd_result_db(checked = c(1L, 0L))
    on.exit(DBI::dbDisconnect(conn), add = TRUE)

    res <- get_xic_traces_from_compounds(conn, NULL)

    expect_equal(nrow(res), 4)
    expect_setequal(unique(res$name), c("Betaine", "Creatine"))
})

test_that("retention times are converted to seconds", {
    conn <- .write_cd_result_db(checked = NULL)
    on.exit(DBI::dbDisconnect(conn), add = TRUE)

    res <- get_xic_traces_from_compounds(conn, 'name == "Betaine"')

    expect_equal(unique(res$rt), 0.86 * 60)
    expect_equal(unique(res$rtmin), 0.80 * 60)
    expect_equal(unique(res$rtmax), 0.92 * 60)
})
