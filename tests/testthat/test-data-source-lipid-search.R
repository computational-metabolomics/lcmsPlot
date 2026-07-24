test_that("LipidSearchSource builds compounds, peaks and metadata", {
    path <- .write_lipid_search_results()
    sample_paths <- file.path(tempdir(), c("ko15.raw", "wt15.raw"))

    res <- LipidSearchSource(path, sample_paths)

    expect_s4_class(res, "LipidSearchSource")
    expect_equal(res@name, "lipid-search")

    expect_equal(nrow(res@metadata), 2)
    expect_true(all(
        c("sample_index", "sample_id", "sample_path", "sample_key",
          "sample_group", "in_results") %in% colnames(res@metadata)
    ))
    expect_true(all(res@metadata$in_results))
    expect_equal(res@metadata$sample_key, c("c-1", "c-2"))
    expect_equal(res@metadata$sample_group, c("c", "c"))

    # Three accepted lipids x two samples.
    expect_equal(nrow(res@compounds), 6)

    # Four peaks were integrated across the accepted lipids.
    expect_equal(nrow(res@peaks), 4)
    expect_true(all(res@peaks$into > 0))
    expect_true(all(
        c("mz", "rt", "rtmin", "rtmax", "into", "maxo", "sample_index",
          "name", "formula", "adduct", "class", "grade") %in%
            colnames(res@peaks)
    ))
})

test_that("LipidSearchSource matches sample paths by basename, ignoring the extension", {
    path <- .write_lipid_search_results(file_names = c("ko15.raw", "wt15.raw"))
    sample_paths <- file.path(tempdir(), c("ko15.mzML", "wt15.mzML"))

    res <- LipidSearchSource(path, sample_paths)

    expect_true(all(res@metadata$in_results))
    expect_equal(res@metadata$sample_key, c("c-1", "c-2"))
})

test_that("LipidSearchSource accepts samples that are absent from the results", {
    path <- .write_lipid_search_results()
    sample_paths <- file.path(
        tempdir(), c("ko15.mzML", "wt15.mzML", "extra.mzML"))

    res <- LipidSearchSource(path, sample_paths)

    expect_equal(nrow(res@metadata), 3)
    expect_equal(res@metadata$in_results, c(TRUE, TRUE, FALSE))
    expect_true(is.na(res@metadata$sample_key[3]))

    extra <- res@compounds[res@compounds$sample_index == 3, ]

    # Every lipid is still extractable in the extra sample...
    expect_equal(nrow(extra), 3)
    expect_false(any(extra$detected))
    expect_true(all(is.na(extra$into)))
    expect_true(all(is.na(extra$maxo)))

    # ...using the consensus m/z and RT of the samples that did detect it.
    # The second lipid was only found in c-2, at m/z 372.3110 and RT 2.25 min.
    acca <- extra[extra$lipid_row == 2, ]
    expect_equal(acca$mz, 372.3110)
    expect_equal(acca$rt, 2.25 * 60)

    # ...and it contributes no peaks.
    expect_false(any(res@peaks$sample_index == 3))
})

test_that("LipidSearchSource reports declared samples with no supplied path", {
    path <- .write_lipid_search_results()

    expect_message(
        res <- LipidSearchSource(path, file.path(tempdir(), "ko15.raw")),
        "wt15.raw"
    )

    expect_equal(nrow(res@metadata), 1)
    expect_equal(res@metadata$sample_key, "c-1")
})

test_that("LipidSearchSource drops rejected lipids unless asked to keep them", {
    path <- .write_lipid_search_results()
    sample_paths <- file.path(tempdir(), c("ko15.raw", "wt15.raw"))

    kept <- LipidSearchSource(path, sample_paths)
    all_rows <- LipidSearchSource(path, sample_paths, keep_rejected = TRUE)

    expect_equal(length(unique(kept@compounds$lipid_row)), 3)
    expect_false(any(grepl("LPC", kept@compounds$lipid_ion)))

    expect_equal(length(unique(all_rows@compounds$lipid_row)), 4)
    expect_true(any(grepl("LPC", all_rows@compounds$lipid_ion)))
})

test_that("LipidSearchSource disambiguates lipid ions reported at several RTs", {
    path <- .write_lipid_search_results()
    sample_paths <- file.path(tempdir(), c("ko15.raw", "wt15.raw"))

    res <- LipidSearchSource(path, sample_paths)
    names <- unique(res@compounds$name)

    expect_equal(length(names), 3)
    expect_true("PC(16:0_18:1)+H" %in% names)
    expect_setequal(
        grep("AcCa", names, value = TRUE),
        c("AcCa(14:0)+H [RT 2.25]", "AcCa(14:0)+H [RT 2.60]")
    )
})

test_that("LipidSearchSource converts RTs to seconds and derives the peak window", {
    path <- .write_lipid_search_results()
    sample_paths <- file.path(tempdir(), c("ko15.raw", "wt15.raw"))

    res <- LipidSearchSource(path, sample_paths)

    pc <- res@compounds[
        res@compounds$lipid_row == 1 & res@compounds$sample_index == 1, ]

    expect_equal(pc$rt, 3.00 * 60)
    expect_equal(pc$rtmin, (3.00 - 0.05) * 60)
    expect_equal(pc$rtmax, (3.00 + 0.06) * 60)
    expect_equal(pc$mz, 760.5855)
    expect_equal(pc$adduct, "+H")
    expect_equal(pc$class, "PC")
    expect_true(pc$detected)
})

test_that("LipidSearchSource falls back to the consensus in samples with no hit", {
    path <- .write_lipid_search_results()
    sample_paths <- file.path(tempdir(), c("ko15.raw", "wt15.raw"))

    res <- LipidSearchSource(path, sample_paths)

    # The second lipid has no observed m/z in c-1, but a per-sample RT.
    acca <- res@compounds[
        res@compounds$lipid_row == 2 & res@compounds$sample_index == 1, ]

    expect_false(acca$detected)
    expect_equal(acca$mz, 372.3110)   # consensus, from c-2
    expect_equal(acca$rt, 2.20 * 60)  # this sample's own RT
    expect_true(is.na(acca$into))
})

test_that("LipidSearchSource errors on a file with neither declarations nor sample columns", {
    path <- tempfile(fileext = ".txt")
    writeLines(c("Rej.\tLipidIon", "0\tPC(16:0)+H"), path)

    expect_error(
        LipidSearchSource(path, file.path(tempdir(), "ko15.raw")),
        "Could not identify any samples"
    )
})

test_that("LipidSearchSource attaches user-supplied metadata", {
    path <- .write_lipid_search_results()
    sample_paths <- file.path(tempdir(), c("ko15.raw", "wt15.raw"))

    res <- LipidSearchSource(
        path,
        sample_paths,
        metadata = data.frame(group = c("KO", "WT"))
    )

    expect_equal(res@metadata$group, c("KO", "WT"))
    expect_error(
        LipidSearchSource(path, sample_paths, metadata = data.frame(group = 1)),
        "one row per sample path"
    )
})

test_that("LipidSearchSource auto-detects 5.2 and keys samples s1, s2", {
    path <- .write_lipid_search5_results()
    sample_paths <- file.path(tempdir(), c("a.mzML", "b.mzML"))

    res <- LipidSearchSource(path, sample_paths)

    expect_s4_class(res, "LipidSearchSource")
    expect_equal(res@version, "5.2")
    expect_equal(nrow(res@metadata), 2)
    expect_equal(res@metadata$sample_key, c("s1", "s2"))
    expect_true(all(res@metadata$in_results))
})

test_that("LipidSearchSource 5.2 maps LipidID, AdductIon and SubClass", {
    path <- .write_lipid_search5_results()
    sample_paths <- file.path(tempdir(), c("a.mzML", "b.mzML"))

    res <- LipidSearchSource(path, sample_paths)

    pc <- res@compounds[
        res@compounds$lipid_row == 1 & res@compounds$sample_index == 1, ]
    expect_equal(pc$name, "PC(16:0_18:1)+H")
    expect_equal(pc$adduct, "M+H")
    expect_equal(pc$sub_class, "diacyl")
    expect_equal(pc$rt, 3.00 * 60)   # ObsRt, minutes -> seconds
    expect_true(pc$detected)
})

test_that("LipidSearchSource 5.2 treats NP as undetected and uses BaseRt as consensus", {
    path <- .write_lipid_search5_results()
    sample_paths <- file.path(tempdir(), c("a.mzML", "b.mzML"))

    res <- LipidSearchSource(path, sample_paths)

    # Row 2 (AcCa) is NP in s1: undetected, no peak, but still extractable via
    # its consensus (observed m/z / RT come from s2, or BaseRt as a fallback).
    acca_s1 <- res@compounds[
        res@compounds$lipid_row == 2 & res@compounds$sample_index == 1, ]
    expect_false(acca_s1$detected)
    expect_true(is.na(acca_s1$into))
    expect_equal(acca_s1$mz, 372.3110)       # consensus obs m/z from s2
    expect_equal(acca_s1$rt, 2.25 * 60)      # consensus obs RT from s2

    expect_true(all(res@peaks$into > 0))
    expect_false(any(is.na(res@peaks$into)))
})

test_that("LipidSearchSource 5.2 rejects Rej == true unless kept", {
    path <- .write_lipid_search5_results()
    sample_paths <- file.path(tempdir(), c("a.mzML", "b.mzML"))

    kept <- LipidSearchSource(path, sample_paths)
    all_rows <- LipidSearchSource(path, sample_paths, keep_rejected = TRUE)

    expect_false(any(grepl("LPC", kept@compounds$lipid_ion)))
    expect_true(any(grepl("LPC", all_rows@compounds$lipid_ion)))
    expect_type(kept@compounds$rej, "logical")
})

test_that("LipidSearchSource 5.2 maps named sample_paths by key regardless of order", {
    path <- .write_lipid_search5_results()
    p1 <- file.path(tempdir(), "a.mzML")
    p2 <- file.path(tempdir(), "b.mzML")

    res <- LipidSearchSource(path, c("s2" = p2, "s1" = p1))

    # The metadata follows sample_paths order (s2 then s1) but the keys are
    # correctly assigned by name.
    expect_equal(res@metadata$sample_key, c("s2", "s1"))
    expect_equal(res@metadata$sample_path, c(p2, p1))
    expect_true(all(res@metadata$in_results))

    # s1's detected PC peak lands on the s1 row (sample_index 2 here).
    s1_index <- res@metadata$sample_index[res@metadata$sample_key == "s1"]
    pc <- res@compounds[
        res@compounds$lipid_row == 1 & res@compounds$sample_index == s1_index, ]
    expect_equal(pc$rt, 3.00 * 60)
})

test_that("LipidSearchSource 5.2 supports a keyed subset plus an unnamed extra", {
    path <- .write_lipid_search5_results()
    p1 <- file.path(tempdir(), "a.mzML")
    extra <- file.path(tempdir(), "extra.mzML")

    expect_message(
        res <- LipidSearchSource(path, c("s1" = p1, extra)),
        "s2"
    )

    expect_equal(res@metadata$sample_key, c("s1", NA))
    expect_equal(res@metadata$in_results, c(TRUE, FALSE))

    # The extra sample still gets every lipid, from the consensus, with no peak.
    extra_rows <- res@compounds[res@compounds$sample_index == 2, ]
    expect_gt(nrow(extra_rows), 0)
    expect_false(any(extra_rows$detected))
    expect_false(any(res@peaks$sample_index == 2))
})

test_that("LipidSearchSource 5.2 maps unnamed sample_paths positionally", {
    path <- .write_lipid_search5_results()
    sample_paths <- file.path(tempdir(), c("first.mzML", "second.mzML"))

    res <- LipidSearchSource(path, sample_paths)

    expect_equal(res@metadata$sample_key, c("s1", "s2"))
    expect_true(all(res@metadata$in_results))
})
