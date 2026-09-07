# Build a minimal Compound Discoverer scripting-node export (node_args.json plus
# the referenced tab-delimited tables) in a temporary directory and return the
# JSON path.
.write_cd_node_export <- function(
    dir = tempfile("cd_node"),
    study_file_ids = c(1, 2),
    file_names = c("ko15.mzML", "wt15.mzML"),
    include_study_files = TRUE,
    checked = NULL
) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)

    write_tab <- function(df, name) {
        path <- file.path(dir, name)
        utils::write.table(
            df, path, sep = "\t", row.names = FALSE, quote = FALSE)
        path
    }

    # Two compounds.
    compounds <- data.frame(
        "Compounds ID" = c(1, 2),
        "Name" = c("Betaine", "Creatine"),
        "Formula" = c("C5 H11 N O2", "C4 H9 N3 O2"),
        check.names = FALSE
    )

    # Compound Discoverer only exports "Checked" once compounds have been
    # checked in the Compounds table, and writes it as the text True/False.
    if (!is.null(checked)) {
        compounds$Checked <- checked
    }

    # Compounds per File: each compound detected in each study file.
    cpf <- data.frame(
        "Compounds per File ID" = c(1, 2, 3, 4),
        "StudyFileID" = c(study_file_ids[1], study_file_ids[2],
                          study_file_ids[1], study_file_ids[2]),
        "Area" = c(5000, 3000, 8000, 6000),
        "Intensity" = c(500, 300, 800, 600),
        check.names = FALSE
    )

    # Features: one molecular ion + one fragment per Compounds-per-File row.
    features <- data.frame(
        "Features ID" = seq_len(8),
        "mz" = c(118.0863, 90.05, 118.0863, 90.05,
                 132.0766, 90.05, 132.0766, 90.05),
        "Ion" = rep(c("[M+H]+1", "[M+H-46]+1"), times = 4),
        "Area" = c(5000, 1000, 3000, 800, 8000, 1500, 6000, 1200),
        "RT [min]" = c(0.86, 0.86, 0.86, 0.86, 0.90, 0.90, 0.90, 0.90),
        "Left RT [min]" = c(0.80, 0.80, 0.80, 0.80, 0.84, 0.84, 0.84, 0.84),
        "Right RT [min]" = c(0.92, 0.92, 0.92, 0.92, 0.96, 0.96, 0.96, 0.96),
        check.names = FALSE
    )

    link_cmp_cpf <- data.frame(
        "Compounds ID" = c(1, 1, 2, 2),
        "Compounds per File ID" = c(1, 2, 3, 4),
        check.names = FALSE
    )

    link_cpf_feat <- data.frame(
        "Compounds per File ID" = c(1, 1, 2, 2, 3, 3, 4, 4),
        "Features ID" = seq_len(8),
        check.names = FALSE
    )

    tables <- list(
        list(DataFile = write_tab(compounds, "compounds.txt")),
        list(DataFile = write_tab(cpf, "cpf.txt")),
        list(DataFile = write_tab(features, "features.txt")),
        list(DataFile = write_tab(link_cmp_cpf, "link_cmp_cpf.txt")),
        list(DataFile = write_tab(link_cpf_feat, "link_cpf_feat.txt"))
    )

    if (include_study_files) {
        study_files <- data.frame(
            "StudyFileID" = study_file_ids,
            "File Name" = file_names,
            check.names = FALSE
        )
        tables <- c(
            tables,
            list(list(DataFile = write_tab(study_files, "study_files.txt")))
        )
    }

    node_args <- list(Tables = tables)
    json_path <- file.path(dir, "node_args.json")
    writeLines(jsonlite::toJSON(node_args, auto_unbox = TRUE), json_path)
    json_path
}

test_that("CompoundDiscovererNodeSource builds peaks and metadata", {
    json_path <- .write_cd_node_export()
    sample_paths <- file.path(tempdir(), c("ko15.mzML", "wt15.mzML"))

    res <- CompoundDiscovererNodeSource(
        node_args = json_path,
        sample_paths = sample_paths
    )

    expect_s4_class(res, "CompoundDiscovererNodeSource")
    expect_equal(res@name, "compound-discoverer-node")

    expect_equal(nrow(res@metadata), 2)
    expect_true(all(
        c("sample_index", "sample_id", "sample_path", "study_file_id") %in%
            colnames(res@metadata)
    ))

    # Two compounds x two samples = four peaks.
    expect_equal(nrow(res@peaks), 4)
    expect_true(all(
        c("mz", "rt", "rtmin", "rtmax", "into", "maxo", "sample_index",
          "name", "formula", "adduct") %in% colnames(res@peaks)
    ))
})

test_that("CompoundDiscovererNodeSource selects the molecular ion m/z", {
    json_path <- .write_cd_node_export()
    sample_paths <- file.path(tempdir(), c("ko15.mzML", "wt15.mzML"))

    res <- CompoundDiscovererNodeSource(json_path, sample_paths)

    betaine <- res@peaks[res@peaks$name == "Betaine", ]
    # [M+H]+1 m/z is selected, not the more/less abundant fragment.
    expect_true(all(abs(betaine$mz - 118.0863) < 1e-4))
    expect_true(all(betaine$adduct == "[M+H]+1"))
})

test_that("CompoundDiscovererNodeSource converts RT minutes to seconds", {
    json_path <- .write_cd_node_export()
    sample_paths <- file.path(tempdir(), c("ko15.mzML", "wt15.mzML"))

    res <- CompoundDiscovererNodeSource(json_path, sample_paths)

    betaine <- res@peaks[res@peaks$name == "Betaine", ][1, ]
    expect_equal(betaine$rt, 0.86 * 60, tolerance = 1e-6)
    expect_equal(betaine$rtmin, 0.80 * 60, tolerance = 1e-6)
    expect_equal(betaine$rtmax, 0.92 * 60, tolerance = 1e-6)
})

test_that("CompoundDiscovererNodeSource errors on unmatched sample basename", {
    json_path <- .write_cd_node_export()
    sample_paths <- file.path(tempdir(), c("other1.mzML", "other2.mzML"))

    expect_error(
        CompoundDiscovererNodeSource(json_path, sample_paths),
        "Could not match study file"
    )
})

test_that("CompoundDiscovererNodeSource matches positionally without a study-file table", {
    json_path <- .write_cd_node_export(include_study_files = FALSE)
    sample_paths <- file.path(tempdir(), c("anyA.mzML", "anyB.mzML"))

    res <- CompoundDiscovererNodeSource(json_path, sample_paths)

    expect_equal(nrow(res@metadata), 2)
    expect_equal(res@metadata$sample_path, sample_paths)
})

# -----------------------------------------------------------------------------
# Compound-centric ("Unknown Compounds") topology: m/z lives on the consolidated
# compound, per-file RT/area live on the instance ("Compounds per File") table,
# ions/features link to the compound, and Name is empty.
# -----------------------------------------------------------------------------
.write_cd_node_unknown_export <- function(
    dir = tempfile("cd_node_unknown"),
    study_file_ids = c("F4", "F5"),
    file_names = c("ko15.mzML", "wt15.mzML"),
    include_features = TRUE
) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)

    write_tab <- function(df, name) {
        path <- file.path(dir, name)
        utils::write.table(
            df, path, sep = "\t", row.names = FALSE, quote = FALSE)
        path
    }

    # Two unknown compounds: empty Name, identified by Formula; m/z is on the
    # consolidated compound. Compound 1 is more abundant than compound 2.
    compounds <- data.frame(
        "Compounds ID" = c(1, 2),
        "Name" = c("", ""),
        "Formula" = c("C4 H7 N", "C3 H5 N O"),
        "m/z" = c(118.0863, 132.0766),
        check.names = FALSE
    )

    # Per-file instances: each compound in each study file, with RT window/area.
    cpf <- data.frame(
        "Compounds per File ID" = c(1, 2, 3, 4),
        "Study File ID" = c(study_file_ids[1], study_file_ids[2],
                            study_file_ids[1], study_file_ids[2]),
        "Apex RT [min]" = c(0.86, 0.86, 0.90, 0.90),
        "Left RT [min]" = c(0.80, 0.80, 0.84, 0.84),
        "Right RT [min]" = c(0.92, 0.92, 0.96, 0.96),
        "Area" = c(8000, 6000, 5000, 3000),
        "Intensity" = c(800, 600, 500, 300),
        check.names = FALSE
    )

    link_cmp_cpf <- data.frame(
        "Compounds ID" = c(1, 1, 2, 2),
        "Compounds per File ID" = c(1, 2, 3, 4),
        check.names = FALSE
    )

    tables <- list(
        list(DataFile = write_tab(compounds, "compounds.txt")),
        list(DataFile = write_tab(cpf, "cpf.txt")),
        list(DataFile = write_tab(link_cmp_cpf, "link_cmp_cpf.txt"))
    )

    if (include_features) {
        # Ion table (compound-level) + Compound<->Features link, for the adduct.
        features <- data.frame(
            "Features ID" = c(1, 2),
            "Ion" = c("[M+H]+1", "[M+H]+1"),
            check.names = FALSE
        )
        link_cmp_feat <- data.frame(
            "Compounds ID" = c(1, 2),
            "Features ID" = c(1, 2),
            check.names = FALSE
        )
        tables <- c(
            tables,
            list(list(DataFile = write_tab(features, "features.txt"))),
            list(list(DataFile = write_tab(link_cmp_feat, "link_cmp_feat.txt")))
        )
    }

    study_files <- data.frame(
        "Study File ID" = study_file_ids,
        "File Name" = file_names,
        check.names = FALSE
    )
    tables <- c(
        tables,
        list(list(DataFile = write_tab(study_files, "study_files.txt")))
    )

    node_args <- list(Tables = tables)
    json_path <- file.path(dir, "node_args.json")
    writeLines(jsonlite::toJSON(node_args, auto_unbox = TRUE), json_path)
    json_path
}

test_that("Compound-centric export: m/z from compound, RT/area from instance", {
    json_path <- .write_cd_node_unknown_export()
    sample_paths <- file.path(tempdir(), c("ko15.mzML", "wt15.mzML"))

    res <- CompoundDiscovererNodeSource(json_path, sample_paths)

    expect_s4_class(res, "CompoundDiscovererNodeSource")
    # Two compounds x two files.
    expect_equal(nrow(res@peaks), 4)

    # m/z comes from the consolidated compound.
    c4 <- res@compounds[res@compounds$formula == "C4 H7 N", ]
    expect_true(all(abs(c4$mz - 118.0863) < 1e-4))

    # RT window comes from the instance table, converted minutes -> seconds.
    expect_equal(c4$rt[1], 0.86 * 60, tolerance = 1e-6)
    expect_equal(c4$rtmin[1], 0.80 * 60, tolerance = 1e-6)
    expect_equal(c4$rtmax[1], 0.92 * 60, tolerance = 1e-6)
    expect_equal(sort(c4$into), c(6000, 8000))
})

test_that("Compound-centric export: names are synthesized, non-empty and unique", {
    json_path <- .write_cd_node_unknown_export()
    sample_paths <- file.path(tempdir(), c("ko15.mzML", "wt15.mzML"))

    res <- CompoundDiscovererNodeSource(json_path, sample_paths)

    expect_true(all(nzchar(res@compounds$name)))
    # One label per compound (shared across its files); two distinct compounds.
    expect_equal(length(unique(res@compounds$name)), 2)
    expect_setequal(unique(res@compounds$name), c("C4 H7 N", "C3 H5 N O"))
})

test_that("Compound-centric export: compound_rank orders by total area", {
    json_path <- .write_cd_node_unknown_export()
    sample_paths <- file.path(tempdir(), c("ko15.mzML", "wt15.mzML"))

    res <- CompoundDiscovererNodeSource(json_path, sample_paths)

    ranks <- unique(res@compounds[, c("formula", "compound_rank")])
    # C4 H7 N total area 14000 > C3 H5 N O total area 8000.
    expect_equal(
        ranks$compound_rank[ranks$formula == "C4 H7 N"], 1)
    expect_equal(
        ranks$compound_rank[ranks$formula == "C3 H5 N O"], 2)

    top1 <- res@compounds[res@compounds$compound_rank <= 1, ]
    expect_setequal(unique(top1$formula), "C4 H7 N")
})

test_that("Compound-centric export: adduct taken from the molecular ion", {
    json_path <- .write_cd_node_unknown_export()
    sample_paths <- file.path(tempdir(), c("ko15.mzML", "wt15.mzML"))

    res <- CompoundDiscovererNodeSource(json_path, sample_paths)

    expect_true(all(res@compounds$adduct == "[M+H]+1"))
})

test_that("Compound-centric export works without a Features table (adduct NA)", {
    json_path <- .write_cd_node_unknown_export(include_features = FALSE)
    sample_paths <- file.path(tempdir(), c("ko15.mzML", "wt15.mzML"))

    res <- CompoundDiscovererNodeSource(json_path, sample_paths)

    expect_equal(nrow(res@peaks), 4)
    expect_true(all(is.na(res@compounds$adduct)))
    # m/z + RT still resolve from compound/instance.
    expect_true(all(res@compounds$mz > 0))
})

# -----------------------------------------------------------------------------
# CD 3.5 "Unknown Compounds" column names as actually exported: "mz",
# "RT in min", "FWHM in min", "Area Ref Ion", "Intensity Max", "Reference Ion".
# The RT window is derived from FWHM (no Left/Right RT columns are exported).
# -----------------------------------------------------------------------------
.write_cd_node_cd35_export <- function(
    dir = tempfile("cd_node_cd35"),
    study_file_ids = c("F4", "F5"),
    file_names = c("ko15.mzML", "wt15.mzML")
) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    write_tab <- function(df, name) {
        path <- file.path(dir, name)
        utils::write.table(
            df, path, sep = "\t", row.names = FALSE, quote = FALSE)
        path
    }

    compounds <- data.frame(
        "Compounds ID" = c(1, 2),
        "Name" = c("", ""),
        "Formula" = c("C4 H7 N", "C3 H5 N O"),
        "mz" = c(70.06509, 72.04437),
        "Reference Ion" = c("[M+H]+1", "[M+H]+1"),
        check.names = FALSE
    )

    cpf <- data.frame(
        "Compounds per File ID" = c(1, 2, 3, 4),
        "Study File ID" = c(study_file_ids[1], study_file_ids[2],
                            study_file_ids[1], study_file_ids[2]),
        "RT in min" = c(7.052, 7.050, 4.498, 4.500),
        "FWHM in min" = c(0.072, 0.072, 0.049, 0.049),
        "Area Ref Ion" = c(8000, 6000, 5000, 3000),
        "Intensity Max" = c(800, 600, 500, 300),
        check.names = FALSE
    )

    link_cmp_cpf <- data.frame(
        "Compounds ID" = c(1, 1, 2, 2),
        "Compounds per File ID" = c(1, 2, 3, 4),
        check.names = FALSE
    )

    study_files <- data.frame(
        "Study File ID" = study_file_ids,
        "File Name" = file_names,
        check.names = FALSE
    )

    tables <- list(
        list(DataFile = write_tab(compounds, "compounds.txt")),
        list(DataFile = write_tab(cpf, "cpf.txt")),
        list(DataFile = write_tab(link_cmp_cpf, "link_cmp_cpf.txt")),
        list(DataFile = write_tab(study_files, "study_files.txt"))
    )

    json_path <- file.path(dir, "node_args.json")
    writeLines(
        jsonlite::toJSON(list(Tables = tables), auto_unbox = TRUE), json_path)
    json_path
}

test_that("CD 3.5 names resolve: mz/rt/into/maxo populated, window from FWHM", {
    json_path <- .write_cd_node_cd35_export()
    sample_paths <- file.path(tempdir(), c("ko15.mzML", "wt15.mzML"))

    res <- CompoundDiscovererNodeSource(json_path, sample_paths)

    expect_equal(nrow(res@peaks), 4)
    # Nothing critical should be NA (this is what caused the extraction crash).
    expect_false(anyNA(res@compounds$mz))
    expect_false(anyNA(res@compounds$rt))
    expect_false(anyNA(res@compounds$rtmin))
    expect_false(anyNA(res@compounds$rtmax))
    expect_false(anyNA(res@compounds$into))
    expect_false(anyNA(res@compounds$maxo))

    c4 <- res@compounds[res@compounds$formula == "C4 H7 N", ][1, ]
    # m/z from the compound; RT from the instance (minutes -> seconds).
    expect_equal(c4$mz, 70.06509, tolerance = 1e-5)
    expect_equal(c4$rt, 7.052 * 60, tolerance = 1e-6)
    # Window derived from FWHM: rt +/- FWHM/2, in seconds.
    expect_equal(c4$rtmin, (7.052 - 0.072 / 2) * 60, tolerance = 1e-6)
    expect_equal(c4$rtmax, (7.052 + 0.072 / 2) * 60, tolerance = 1e-6)
    expect_equal(c4$into, 8000)
    expect_equal(c4$maxo, 800)
})

test_that("CD 3.5 names: adduct comes from the compound Reference Ion", {
    json_path <- .write_cd_node_cd35_export()
    sample_paths <- file.path(tempdir(), c("ko15.mzML", "wt15.mzML"))

    res <- CompoundDiscovererNodeSource(json_path, sample_paths)

    expect_true(all(res@compounds$adduct == "[M+H]+1"))
})

test_that("compound_id (Compounds ID) is retained on the compounds table", {
    # Compound-centric fixture.
    json_path <- .write_cd_node_unknown_export()
    sample_paths <- file.path(tempdir(), c("ko15.mzML", "wt15.mzML"))
    res <- CompoundDiscovererNodeSource(json_path, sample_paths)

    expect_true("compound_id" %in% colnames(res@compounds))
    # Two compounds, ids 1 and 2; each label maps to exactly one id.
    id_by_name <- unique(res@compounds[, c("name", "compound_id")])
    expect_setequal(id_by_name$compound_id, c(1, 2))
    expect_equal(nrow(id_by_name), 2)

    # Per-file-feature fixture retains it too.
    json2 <- .write_cd_node_export()
    res2 <- CompoundDiscovererNodeSource(
        json2, file.path(tempdir(), c("ko15.mzML", "wt15.mzML")))
    expect_true("compound_id" %in% colnames(res2@compounds))
    expect_setequal(unique(res2@compounds$compound_id), c(1, 2))
})

# -----------------------------------------------------------------------------
# node_response.json writer: adds a per-compound filename column to Compounds.
# -----------------------------------------------------------------------------
test_that(".cd_node_write_plot_column adds a filename column + response JSON", {
    # Build a minimal export with an ExpectedResponsePath in a temp dir.
    dir <- tempfile("cd_resp")
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)

    compounds <- data.frame(
        "Compounds ID" = c(1, 2, 3),
        "Name" = c("", "", ""),
        "Formula" = c("C4 H7 N", "C3 H5 N O", "C2 H6 O"),
        check.names = FALSE
    )
    cmp_file <- file.path(dir, "ConsolidatedUnknownCompoundItem.txt")
    utils::write.table(compounds, cmp_file, sep = "\t", row.names = FALSE)

    # A second table CD would re-validate. Its column name has a trailing space
    # that (as in real CD exports) must NOT end up in the response.
    other_file <- file.path(dir, "WorkflowInputFile.txt")
    utils::write.table(
        data.frame("Input Files " = c(1, 2), check.names = FALSE),
        other_file, sep = "\t", row.names = FALSE)

    response_path <- file.path(dir, "node_response.json")
    node_args <- list(
        ExpectedResponsePath = response_path,
        NodeParameters = setNames(list(), character(0)),
        Tables = list(
            list(
                TableName = "Compounds",
                DataFile = cmp_file,
                DataFormat = "CSV",
                Options = setNames(list(), character(0)),
                ColumnDescriptions = list(
                    list(ColumnName = "Compounds ID", ID = "ID",
                         DataType = "Int",
                         Options = setNames(list(), character(0))),
                    list(ColumnName = "Name", ID = "", DataType = "String",
                         Options = setNames(list(), character(0)))
                )
            ),
            list(
                TableName = "Input Files",
                DataFile = other_file,
                DataFormat = "CSV",
                Options = setNames(list(), character(0)),
                ColumnDescriptions = list(
                    list(ColumnName = "Input Files ", ID = "", DataType = "Int",
                         Options = setNames(list(), character(0)))
                )
            )
        )
    )
    node_args_path <- file.path(dir, "node_args.json")
    writeLines(
        jsonlite::toJSON(node_args, auto_unbox = TRUE, pretty = TRUE),
        node_args_path)

    # Map ids 1 and 3 to plot files (id 2 left unmapped -> blank).
    id_to_value <- data.frame(
        compound_id = c(1, 3),
        value = c("C:/plots/1_C4H7N.png", "C:/plots/3_C2H6O.png"),
        stringsAsFactors = FALSE
    )

    out <- lcmsPlot:::.cd_node_write_plot_column(
        node_args = node_args_path,
        id_to_value = id_to_value,
        column_name = "Plot",
        renderer = "EB29D794-4F2E-4785-8B80-A24D8C0FB3E4",
        position_after = "Name"
    )
    expect_equal(out, response_path)
    expect_true(file.exists(response_path))

    # The response must parse and contain ONLY the modified Compounds table.
    resp <- jsonlite::fromJSON(response_path, simplifyVector = FALSE)
    expect_equal(length(resp$Tables), 1)
    cmp <- resp$Tables[[1]]
    expect_equal(cmp$TableName, "Compounds")
    expect_match(cmp$DataFile, "\\.out\\.txt$")
    expect_true(file.exists(cmp$DataFile))

    # The new column description carries the renderer + type.
    last_col <- cmp$ColumnDescriptions[[length(cmp$ColumnDescriptions)]]
    expect_equal(last_col$ColumnName, "Plot")
    expect_equal(last_col$DataType, "String")
    expect_equal(last_col$Options$SpecialCellRenderer,
                 "EB29D794-4F2E-4785-8B80-A24D8C0FB3E4")
    expect_equal(last_col$Options$PositionAfter, "Name")

    # The out.txt has the Plot column filled by Compounds ID (blank for id 2).
    out_df <- utils::read.table(
        cmp$DataFile, header = TRUE, sep = "\t", check.names = FALSE,
        stringsAsFactors = FALSE, colClasses = "character")
    expect_true("Plot" %in% colnames(out_df))
    expect_equal(out_df$Plot[out_df$`Compounds ID` == "1"],
                 "C:/plots/1_C4H7N.png")
    expect_equal(out_df$Plot[out_df$`Compounds ID` == "3"],
                 "C:/plots/3_C2H6O.png")
    expect_equal(out_df$Plot[out_df$`Compounds ID` == "2"], "")

    # Empty Options must serialize as {} not [].
    raw <- paste(readLines(response_path), collapse = "\n")
    expect_false(grepl('"Options"\\s*:\\s*\\[', raw))
    expect_false(grepl('"NodeParameters"\\s*:\\s*\\[', raw))
})

test_that(".cd_node_write_plot_column returns NULL without ExpectedResponsePath", {
    dir <- tempfile("cd_resp_none")
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    cmp_file <- file.path(dir, "c.txt")
    utils::write.table(
        data.frame("Compounds ID" = 1, "Name" = "x", check.names = FALSE),
        cmp_file, sep = "\t", row.names = FALSE)
    node_args <- list(Tables = list(list(
        TableName = "Compounds", DataFile = cmp_file,
        ColumnDescriptions = list())))
    p <- file.path(dir, "node_args.json")
    writeLines(jsonlite::toJSON(node_args, auto_unbox = TRUE), p)

    out <- lcmsPlot:::.cd_node_write_plot_column(
        p, data.frame(compound_id = 1, value = "a.png"))
    expect_null(out)
})

test_that(".cd_node_parse_checked() handles the encodings CD may export", {
    expect_equal(
        .cd_node_parse_checked(c("True", "False")), c(TRUE, FALSE))
    expect_equal(
        .cd_node_parse_checked(c("TRUE", "false", "")), c(TRUE, FALSE, FALSE))
    expect_equal(.cd_node_parse_checked(c("1", "0")), c(TRUE, FALSE))
    expect_equal(.cd_node_parse_checked(c("Yes", "No")), c(TRUE, FALSE))
    expect_equal(.cd_node_parse_checked(c(1L, 0L)), c(TRUE, FALSE))
    expect_equal(.cd_node_parse_checked(c(TRUE, FALSE)), c(TRUE, FALSE))
    # Anything unrecognised stays NA rather than being read as unchecked.
    expect_true(is.na(.cd_node_parse_checked("???")))
})

test_that("the node source exposes the exported Checked column", {
    json_path <- .write_cd_node_export(checked = c("True", "False"))
    sample_paths <- file.path(tempdir(), c("ko15.mzML", "wt15.mzML"))

    res <- CompoundDiscovererNodeSource(json_path, sample_paths)

    expect_true("checked" %in% colnames(res@compounds))
    expect_type(res@compounds$checked, "logical")

    # Betaine is compound 1 (checked), Creatine compound 2 (unchecked); each is
    # detected in both study files.
    expect_true(all(res@compounds$checked[res@compounds$name == "Betaine"]))
    expect_false(any(res@compounds$checked[res@compounds$name == "Creatine"]))

    # The check state is compound metadata, not peak data.
    expect_false("checked" %in% colnames(res@peaks))
})

test_that("the node source omits checked when the export has no such column", {
    json_path <- .write_cd_node_export()
    sample_paths <- file.path(tempdir(), c("ko15.mzML", "wt15.mzML"))

    res <- CompoundDiscovererNodeSource(json_path, sample_paths)

    expect_false("checked" %in% colnames(res@compounds))
})

test_that("compounds can be filtered on checked in the node source", {
    json_path <- .write_cd_node_export(checked = c("True", "False"))
    sample_paths <- file.path(tempdir(), c("ko15.mzML", "wt15.mzML"))

    res <- CompoundDiscovererNodeSource(json_path, sample_paths)

    kept <- dplyr::filter(res@compounds, .data$checked)
    expect_equal(unique(kept$name), "Betaine")

    dropped <- dplyr::filter(res@compounds, !.data$checked)
    expect_equal(unique(dropped$name), "Creatine")
})
