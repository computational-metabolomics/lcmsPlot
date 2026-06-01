test_that(".mz_window returns symmetric ppm window around mz", {
    w <- .mz_window(500, 10)
    expect_equal(length(w), 2)
    expect_equal(w[1], 500 - 500 * 10 / 1e6)
    expect_equal(w[2], 500 + 500 * 10 / 1e6)
    expect_lt(w[1], w[2])
})

test_that(".mz_window scales linearly with ppm", {
    expect_equal(diff(.mz_window(100, 20)), 2 * 100 * 20 / 1e6)
    expect_equal(diff(.mz_window(1000, 5)),  2 * 1000 * 5 / 1e6)
})

test_that(".eic_features returns NA RT bounds when rt is NULL", {
    m <- .eic_features(mz = 335, ppm = 20, rt = NULL)
    expect_true(is.matrix(m))
    expect_equal(nrow(m), 1)
    expect_equal(colnames(m), c("mzmin", "mzmax", "rtmin", "rtmax"))
    expect_true(is.na(m[1, "rtmin"]))
    expect_true(is.na(m[1, "rtmax"]))
})

test_that(".eic_features fills RT window when rt and rt_tol supplied", {
    m <- .eic_features(mz = 335, ppm = 20, rt = 2800, rt_tol = 100)
    expect_equal(unname(m[1, "rtmin"]), 2700)
    expect_equal(unname(m[1, "rtmax"]), 2900)
})

test_that(".eic_features treats NA rt the same as NULL", {
    m <- .eic_features(mz = 335, ppm = 20, rt = NA_real_, rt_tol = 100)
    expect_true(is.na(m[1, "rtmin"]))
    expect_true(is.na(m[1, "rtmax"]))
})

test_that(".uploaded_files_to_paths returns empty for NULL or empty upload", {
    expect_identical(
        .uploaded_files_to_paths(NULL, tempfile()),
        character(0))
    empty <- data.frame(name = character(), datapath = character())
    expect_identical(
        .uploaded_files_to_paths(empty, tempfile()),
        character(0))
})

test_that(".load_dataset returns mzML/CDF/raw paths unchanged", {
    paths <- c("/tmp/a.mzML", "/tmp/b.CDF", "/tmp/c.raw")
    expect_identical(.load_dataset(paths), paths)

    expect_identical(
        .load_dataset(c("/tmp/x.MZML", "/tmp/y.cdf")),
        c("/tmp/x.MZML", "/tmp/y.cdf"))
})

test_that(".load_dataset returns a .cdResult path unchanged", {
    expect_identical(.load_dataset("/tmp/foo.cdResult"),
                     "/tmp/foo.cdResult")
})

test_that(".load_dataset readRDS-loads a .rds upload", {
    f <- tempfile(fileext = ".rds")
    saveRDS(list(x = 1, y = "hello"), f)
    out <- .load_dataset(f)
    expect_equal(out, list(x = 1, y = "hello"))
    unlink(f)
})

test_that(".load_dataset extracts a supported object from an .RData", {
    skip_if_not_installed("xcms")
    fake <- structure(list(), class = "XCMSnExp")
    f <- tempfile(fileext = ".RData")
    save(fake, file = f)
    out <- .load_dataset(f)
    expect_s3_class(out, "XCMSnExp")
    unlink(f)
})

test_that(".load_dataset errors when .RData has no supported object", {
    junk <- list(a = 1, b = 2)
    f <- tempfile(fileext = ".RData")
    save(junk, file = f)
    expect_error(.load_dataset(f), "does not contain")
    unlink(f)
})

test_that(".load_dataset rejects mixed raw + rds uploads", {
    expect_error(
        .load_dataset(c("/tmp/a.mzML", "/tmp/b.rds")),
        "Only one file")
})

test_that(".load_dataset rejects unknown extensions", {
    expect_error(.load_dataset("/tmp/foo.xls"), "Unsupported file type")
})

test_that(".load_dataset rejects empty input", {
    expect_error(.load_dataset(character(0)), "No files")
})

test_that(".uploaded_files_to_paths copies files under their original names", {
    tmp_src <- tempfile(fileext = ".dat")
    writeLines("hello", tmp_src)
    upload <- data.frame(
        name     = "original.mzML",
        datapath = tmp_src,
        stringsAsFactors = FALSE)

    dest_dir <- tempfile()
    out <- .uploaded_files_to_paths(upload, dest_dir)

    expect_length(out, 1)
    expect_true(file.exists(out))
    expect_equal(basename(out), "original.mzML")
    expect_equal(readLines(out), "hello")

    unlink(tmp_src)
    unlink(dest_dir, recursive = TRUE)
})
