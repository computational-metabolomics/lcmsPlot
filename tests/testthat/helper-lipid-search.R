# Build a minimal LipidSearch 4.2 result file in a temporary location and
# return its path. `lipids` is a data frame of lipid-level fields; `values` is a
# named list of per-sample matrices (one column per sample, one row per lipid).
.write_lipid_search_results <- function(
    path = tempfile(fileext = ".txt"),
    file_names = c("ko15.raw", "wt15.raw"),
    lipids = NULL,
    values = NULL
) {
    sample_keys <- paste0("c-", seq_along(file_names))

    if (is.null(lipids)) {
        lipids <- data.frame(
            "Rej." = c(0, 0, 0, 1),
            "LipidIon" = c(
                "PC(16:0_18:1)+H", "AcCa(14:0)+H", "AcCa(14:0)+H",
                "LPC(18:0)+H"
            ),
            "LipidGroup" = c(
                "PC(16:0_18:1)", "AcCa(14:0)", "AcCa(14:0)", "LPC(18:0)"
            ),
            "Class" = c("PC", "AcCa", "AcCa", "LPC"),
            "FattyAcid" = c("(16:0/18:1)", "(14:0)", "(14:0)", "(18:0)"),
            "CalcMz" = c(760.5851, 372.3108, 372.3108, 524.3711),
            "IonFormula" = c(
                "C42 H83 O8 N1 P1", "C21 H42 O4 N1", "C21 H42 O4 N1",
                "C26 H55 O7 N1 P1"
            ),
            check.names = FALSE
        )
    }

    if (is.null(values)) {
        # Row 1 found in both samples; row 2 only in c-2; row 3 only in c-1;
        # row 4 (rejected) found nowhere.
        values <- list(
            "Area" = rbind(c(5000, 3000), c(0, 8000), c(1200, 0), c(0, 0)),
            "Height" = rbind(c(500, 300), c(0, 800), c(120, 0), c(0, 0)),
            "Rt" = rbind(
                c(3.00, 3.02), c(2.20, 2.25), c(2.60, 2.62), c(4.10, 4.12)),
            "ObsMz" = rbind(
                c("760.5855", "760.5849"), c("", "372.3110"),
                c("372.3105", ""), c("", "")),
            "Hwhm(L)" = rbind(
                c(0.05, 0.04), c(0, 0.02), c(0.03, 0), c(0, 0)),
            "Hwhm(R)" = rbind(
                c(0.06, 0.05), c(0, 0.03), c(0.03, 0), c(0, 0)),
            "Grade" = rbind(
                c("A", "B"), c("", "A"), c("C", ""), c("", "")),
            "mScore" = rbind(
                c("8", "7"), c("", "9"), c("6", ""), c("", "")),
            "S/N" = rbind(
                c("20", "15"), c("", "30"), c("12", ""), c("", ""))
        )
    }

    # Per-sample columns are grouped by base name, e.g. Area[c-1], Area[c-2],
    # Height[c-1], ... - as in a real export.
    per_sample <- do.call(cbind, lapply(names(values), function(base) {
        mat <- as.data.frame(values[[base]], stringsAsFactors = FALSE)
        colnames(mat) <- paste0(base, "[", sample_keys, "]")
        mat
    }))

    table <- cbind(lipids, per_sample)

    header <- c(
        paste0("#[", sample_keys, "]:", file_names),
        "#normalize base:",
        "#control group:c",
        "#toprank filter:true",
        "#mScoreThreshold:5.0",
        ""
    )

    # A real export ends every row - including the header - with a tab.
    body <- c(
        paste0(paste(colnames(table), collapse = "\t"), "\t"),
        apply(table, 1, function(row) {
            paste0(paste(as.character(row), collapse = "\t"), "\t")
        })
    )

    writeLines(c(header, body), path)
    path
}

# Build a minimal LipidSearch 5.2 result file (no #[key]:file declarations; the
# table starts on line 1, samples are keyed s1, s2, ... from OrgMeanArea[...],
# and per-injection data is keyed s1-1, s2-1, ...). Returns its path.
.write_lipid_search5_results <- function(
    path = tempfile(fileext = ".txt"),
    lipids = NULL,
    values = NULL
) {
    sample_keys <- c("s1", "s2")
    injection_keys <- c("s1-1", "s2-1")

    if (is.null(lipids)) {
        # Row 3 is rejected; the same LipidID appears at two RTs (rows 1 & 4).
        lipids <- data.frame(
            "Rej" = c("false", "false", "true", "false"),
            "ID" = c(1, 2, 3, 4),
            "LipidID" = c(
                "PC(16:0_18:1)+H", "AcCa(14:0)+H", "LPC(18:0)+H",
                "AcCa(14:0)+H"),
            "LipidGroup" = c(
                "PC(16:0_18:1)", "AcCa(14:0)", "LPC(18:0)", "AcCa(14:0)"),
            "Charge" = c(1, 1, 1, 1),
            "CalcMz" = c(760.5851, 372.3108, 524.3711, 372.3108),
            "BaseRt" = c(3.00, 2.20, 4.10, 2.60),
            "Class" = c("PC", "AcCa", "LPC", "AcCa"),
            "SubClass" = c("diacyl", "carnitines", "monoacyl", "carnitines"),
            "AdductIon" = c("M+H", "M+H", "M+H", "M+H"),
            "IonFormula" = c(
                "C42 H83 O8 N1 P1", "C21 H42 O4 N1", "C26 H55 O7 N1 P1",
                "C21 H42 O4 N1"),
            check.names = FALSE
        )
    }

    if (is.null(values)) {
        # "NP" marks a lipid not present in a sample. Row 1 detected in both;
        # row 2 only in s2; row 4 only in s1; row 3 (rejected) nowhere.
        NP <- "NP"
        values <- list(
            "ObsMz" = rbind(
                c("760.5855", "760.5849"), c(NP, "372.3110"),
                c(NP, NP), c("372.3105", NP)),
            "ObsRt" = rbind(
                c("3.00", "3.02"), c(NP, "2.25"),
                c(NP, NP), c("2.60", NP)),
            "Grade" = rbind(
                c("A", "B"), c(NP, "A"), c(NP, NP), c("C", NP)),
            "Area" = rbind(
                c("5000", "3000"), c(NP, "8000"), c(NP, NP), c("1200", NP)),
            "Height" = rbind(
                c("500", "300"), c(NP, "800"), c(NP, NP), c("120", NP)),
            "IDScore" = rbind(
                c("8", "7"), c(NP, "9"), c(NP, NP), c("6", NP))
        )
    }

    # Group-level OrgMeanArea columns define the sample keys s1, s2.
    org <- as.data.frame(
        matrix(0, nrow = nrow(lipids), ncol = length(sample_keys)),
        stringsAsFactors = FALSE)
    colnames(org) <- paste0("OrgMeanArea[", sample_keys, "]")

    # Per-injection columns keyed s1-1, s2-1.
    per_injection <- do.call(cbind, lapply(names(values), function(base) {
        mat <- as.data.frame(values[[base]], stringsAsFactors = FALSE)
        colnames(mat) <- paste0(base, "[", injection_keys, "]")
        mat
    }))

    table <- cbind(lipids, org, per_injection)

    writeLines(
        c(
            paste(colnames(table), collapse = "\t"),
            apply(table, 1, function(row) paste(row, collapse = "\t"))
        ),
        path
    )
    path
}
