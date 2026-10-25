test_that("MzTabMSource reads the trimmed HUPO-PSI MTBLS263 example", {
    # arrange
    path <- test_path("fixtures", "mztab", "MTBLS263-trimmed.mztab")
    sample_paths <- file.path(
        tempdir(),
        c(paste0("3injections_inj", 1:3, "_POS.mzML"),
          paste0("3samples_sampl", 1:3, "_POS.mzML")))

    # act
    res <- MzTabMSource(path, sample_paths)

    # assert
    expect_s4_class(res, "MzTabMSource")
    expect_s4_class(res, "ExternalDataSource")
    expect_equal(res@name, "mztab-m")

    expect_true(all(
        c("mz", "rt", "rtmin", "rtmax", "into", "maxo", "sample_index",
          "feature_id", "adduct_ion", "charge") %in% colnames(res@peaks)))

    # 4 features, each with an abundance in all 6 assays.
    expect_equal(nrow(res@peaks), 24)
    expect_setequal(res@peaks$sample_index, 1:6)
    expect_true(all(is.na(res@peaks$maxo)))

    proline <- res@peaks |> dplyr::filter(feature_id == "7", sample_index == 1)
    expect_equal(proline$mz, 116.0699)
    expect_equal(proline$rt, 607.75)
    expect_equal(proline$rtmin, 587.39)
    expect_equal(proline$rtmax, 648.11)
    expect_equal(proline$into, 59134552.3)
    expect_equal(proline$name, "Proline")
    expect_equal(proline$formula, "C5H9NO2")
    expect_equal(proline$reliability, 2)
})

test_that("MzTabMSource attaches one SML row to every feature it references", {
    # arrange
    path <- test_path("fixtures", "mztab", "MTBLS263-trimmed.mztab")
    sample_paths <- file.path(
        tempdir(),
        c(paste0("3injections_inj", 1:3, "_POS.mzML"),
          paste0("3samples_sampl", 1:3, "_POS.mzML")))

    # act
    res <- MzTabMSource(path, sample_paths)
    features <- res@compounds |> dplyr::distinct(feature_id, .keep_all = TRUE)

    # assert - SMF 6 and 937 are two adducts of creatinine ("6 | 937").
    creatinine <- features |> dplyr::filter(feature_id %in% c("6", "937"))
    expect_equal(creatinine$name, c("Creatinine", "Creatinine"))
    expect_equal(creatinine$sml_id, c("469", "469"))
    expect_false(any(creatinine$ambiguous))

    # SMF 13 is referenced by both leucine and isoleucine; the first is kept.
    leucine <- features |> dplyr::filter(feature_id == "13")
    expect_equal(leucine$name, "Leucine")
    expect_true(leucine$ambiguous)
})

test_that("MzTabMSource matches URL-encoded locations by basename", {
    # arrange
    path <- .write_mztab(locations = c(
        "file:///D:/Data%20Sets/run%20one.mzML",
        "file:///C:\\data\\run two.mzML"))
    sample_paths <- file.path(tempdir(), c("run two.mzML", "run one.mzML"))

    # act
    res <- MzTabMSource(path, sample_paths)

    # assert - assay 1 is "run one", the second sample path.
    feature_1 <- res@peaks |> dplyr::filter(feature_id == "1")
    expect_equal(
        feature_1$into[order(feature_1$sample_index)], c(2000, 1000))
})

test_that("MzTabMSource ignores the extension when no basename matches", {
    # arrange
    path <- .write_mztab(locations = c("file:///ko15.raw", "file:///wt15.raw"))
    sample_paths <- file.path(tempdir(), c("ko15.mzML", "wt15.mzML"))

    # act
    res <- MzTabMSource(path, sample_paths)

    # assert
    expect_setequal(res@peaks$sample_index, c(1, 2))
})

test_that("MzTabMSource reports the runs that match no sample path", {
    # arrange
    path <- .write_mztab()
    sample_paths <- file.path(tempdir(), "ko15.mzML")

    # act, assert
    expect_error(
        MzTabMSource(path, sample_paths),
        "Could not match these mzTab-M ms_run locations to sample_paths: wt15.mzML"
    )
})

test_that("MzTabMSource keeps samples without an abundance out of the peaks", {
    # arrange
    path <- .write_mztab()
    sample_paths <- file.path(tempdir(), c("ko15.mzML", "wt15.mzML"))

    # act
    res <- MzTabMSource(path, sample_paths)

    # assert - feature 2 has no abundance in assay 2, feature 3 none in assay 1.
    expect_equal(nrow(res@peaks), 4)
    expect_false(any(
        res@peaks$feature_id == "2" & res@peaks$sample_index == 2))

    # Every feature can still be extracted from every sample.
    expect_equal(nrow(res@compounds), 6)
    undetected <- res@compounds |> dplyr::filter(!detected)
    expect_equal(nrow(undetected), 2)
    expect_true(all(is.na(undetected$into)))
})

test_that("MzTabMSource falls back to the apex when RT bounds are null", {
    # arrange
    path <- .write_mztab()
    sample_paths <- file.path(tempdir(), c("ko15.mzML", "wt15.mzML"))

    # act
    res <- MzTabMSource(path, sample_paths)

    # assert
    feature_3 <- res@peaks |> dplyr::filter(feature_id == "3")
    expect_equal(feature_3$rtmin, 3450)
    expect_equal(feature_3$rtmax, 3450)
})

test_that("MzTabMSource leaves features without an SML row unannotated", {
    # arrange
    path <- .write_mztab(sml = .mztab_default_sml())
    sample_paths <- file.path(tempdir(), c("ko15.mzML", "wt15.mzML"))

    # act
    res <- MzTabMSource(path, sample_paths)
    features <- res@compounds |> dplyr::distinct(feature_id, .keep_all = TRUE)

    # assert
    expect_equal(features$name, c("Compound A", "Compound A", NA))
    expect_equal(features$ambiguous, c(FALSE, FALSE, FALSE))
    expect_equal(features$adduct_ion, c("[M+H]1+", "[M+Na]1+", NA))
})

test_that("MzTabMSource reads files without a summary section", {
    # arrange
    path <- .write_mztab()
    sample_paths <- file.path(tempdir(), c("ko15.mzML", "wt15.mzML"))

    # act
    res <- MzTabMSource(path, sample_paths)

    # assert
    expect_true(all(is.na(res@compounds$name)))
})

test_that("MzTabMSource fails without a feature section", {
    # arrange
    path <- tempfile(fileext = ".mzTab")
    writeLines(c(
        "MTD\tms_run[1]-location\tfile:///ko15.mzML",
        "MTD\tassay[1]-ms_run_ref\tms_run[1]"), path)

    # act, assert
    expect_error(
        MzTabMSource(path, file.path(tempdir(), "ko15.mzML")),
        "no small molecule feature"
    )
})

test_that("MzTabMSource adds a sample_id to the metadata", {
    # arrange
    path <- .write_mztab()
    sample_paths <- file.path(tempdir(), c("ko15.mzML", "wt15.mzML"))

    # act
    res <- MzTabMSource(path, sample_paths)

    # assert
    expect_equal(res@metadata$sample_id, c("ko15", "wt15"))
    expect_equal(
        get_metadata(res, "sample_id", NULL)$sample_index, c(1, 2))
    expect_identical(get_detected_peaks(res), res@peaks)
})

test_that("lp_mztab requires an mzTab-M source", {
    # arrange
    p <- lcmsPlot(get_test_sample_paths(1))

    # act, assert
    expect_error(p + lp_mztab(), "not an mzTab-M source")
})
