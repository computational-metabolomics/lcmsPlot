# Build a minimal mzTab-M 2.0 file in a temporary location and return its path.
# `locations` gives one ms_run location (and assay) per entry; `sml` and `smf`
# are data frames written as the SML and SMF sections, with NA written as the
# literal `null`. A NULL `sml` omits the summary section.
.write_mztab <- function(
    path = tempfile(fileext = ".mzTab"),
    locations = c("file:///data/ko15.mzML", "file:///data/wt15.mzML"),
    sml = NULL,
    smf = NULL
) {
    if (is.null(smf)) {
        smf <- data.frame(
            SMF_ID = c(1, 2, 3),
            SME_ID_REFS = c(1, 2, 3),
            adduct_ion = c("[M+H]1+", "[M+Na]1+", NA),
            exp_mass_to_charge = c(335.0, 357.0, 344.1),
            charge = c(1, 1, 1),
            retention_time_in_seconds = c(2800, 2800, 3450),
            retention_time_in_seconds_start = c(2780, 2780, NA),
            retention_time_in_seconds_end = c(2820, 2820, NA),
            "abundance_assay[1]" = c(1000, 500, NA),
            "abundance_assay[2]" = c(2000, NA, 300),
            check.names = FALSE
        )
    }

    as_section <- function(df, header_prefix, row_prefix) {
        values <- as.matrix(format(df, trim = TRUE, scientific = FALSE))
        values[is.na(df)] <- "null"
        c(
            paste(c(header_prefix, colnames(df)), collapse = "\t"),
            # Real exports often end rows with stray tabs.
            paste0(apply(values, 1, function(row) {
                paste(c(row_prefix, row), collapse = "\t")
            }), "\t\t"),
            ""
        )
    }

    n <- seq_along(locations)
    mtd <- c(
        "MTD\tmzTab-version\t2.0.0-M",
        paste0("MTD\tms_run[", n, "]-location\t", locations),
        paste0("MTD\tassay[", n, "]-ms_run_ref\tms_run[", n, "]"),
        ""
    )

    writeLines(
        c(
            mtd,
            "COM\tSummary rows",
            if (!is.null(sml)) as_section(sml, "SMH", "SML"),
            as_section(smf, "SFH", "SMF")
        ),
        path
    )
    path
}

# SML rows for the default SMF section: one compound seen as two adducts
# (features 1 and 2), and feature 3 left unannotated.
.mztab_default_sml <- function() {
    data.frame(
        SML_ID = 1,
        SMF_ID_REFS = "1 | 2",
        database_identifier = "CHEBI:1",
        chemical_formula = "C10H10",
        chemical_name = "Compound A",
        reliability = 2
    )
}
