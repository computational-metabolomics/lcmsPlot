test_that("MZmineFeatureListsSource creates ExternalDataSource (MZmine 2)", {
    # arrange
    tmp <- tempdir()

    sample_paths <- file.path(tmp, c("sample1.mzML"))
    file.create(sample_paths)

    feature_path <- file.path(tmp, "features.csv")

    write.csv(
        data.frame(
            "row m/z" = 100,
            "sample1.mzML Feature status" = "DETECTED",
            "sample1.mzML Feature m/z" = 100,
            "sample1.mzML Feature RT" = 5,
            "sample1.mzML Peak area" = 1000,
            "sample1.mzML Peak height" = 500,
            "sample1.mzML Feature m/z min" = 99.9,
            "sample1.mzML Feature m/z max" = 100.1,
            "sample1.mzML Feature RT start" = 4.8,
            "sample1.mzML Feature RT end" = 5.2,
            check.names = FALSE
        ),
        feature_path,
        row.names = FALSE
    )

    # act
    res <- MZmineFeatureListsSource(
        feature_lists_paths = feature_path,
        sample_paths = sample_paths
    )

    # assert
    expect_s4_class(res, "ExternalDataSource")
    expect_equal(res@name, "mzmine")

    expect_true(is.data.frame(res@metadata))
    expect_true(is.data.frame(res@peaks))

    expect_equal(nrow(res@peaks), 1)
    expect_true(all(
        c("mz", "rt", "into", "maxo",
          "mzmin", "mzmax", "rtmin", "rtmax",
          "sample_index") %in% colnames(res@peaks)
    ))

    expect_equal(res@peaks$sample_index, 1)
})

test_that("MZmineFeatureListsSource works for MZmine >=3 format", {
    # arrange
    tmp <- tempdir()

    sample_paths <- file.path(tmp, "sampleA.mzML")
    file.create(sample_paths)

    feature_path <- file.path(tmp, "features.tsv")

    write.table(
        data.frame(
            mz = 150,
            rt = 10,
            area = 2000,
            height = 800,
            "mz_range:min" = 149.9,
            "mz_range:max" = 150.1,
            "rt_range:min" = 9.8,
            "rt_range:max" = 10.2,
            "datafile:sampleA.mzML:feature_state" = "DETECTED",
            check.names = FALSE
        ),
        feature_path,
        sep = "\t",
        row.names = FALSE
    )

    # act
    res <- MZmineFeatureListsSource(
        feature_lists_paths = feature_path,
        sample_paths = sample_paths
    )

    # assert
    expect_equal(nrow(res@peaks), 1)
    expect_equal(res@peaks$sample_index, 1)
})

test_that("MZmineFeatureListsSource fails when feature list filename does not match sample_paths", {
    # arrange
    tmp <- tempdir()

    sample_paths <- file.path(tmp, "sample1.mzML")
    file.create(sample_paths)

    feature_path <- file.path(tmp, "features.csv")

    write.csv(
        data.frame(
            "row m/z" = 100,
            "wrongname.mzML Feature status" = "DETECTED",
            "wrongname.mzML Feature m/z" = 100,
            "wrongname.mzML Feature RT" = 5,
            "wrongname.mzML Peak area" = 1000,
            "wrongname.mzML Peak height" = 500,
            "wrongname.mzML Feature m/z min" = 99.9,
            "wrongname.mzML Feature m/z max" = 100.1,
            "wrongname.mzML Feature RT start" = 4.8,
            "wrongname.mzML Feature RT end" = 5.2,
            check.names = FALSE
        ),
        feature_path,
        row.names = FALSE
    )

    # act, assert
    expect_error(
        MZmineFeatureListsSource(feature_path, sample_paths),
        "Could not match feature list file to sample_paths"
    )
})

test_that("MZmineFeatureListsSource fails when required columns are missing", {
    # arrange
    tmp <- tempdir()

    sample_paths <- file.path(tmp, "sample1.mzML")
    file.create(sample_paths)

    feature_path <- file.path(tmp, "features.csv")

    write.csv(
        data.frame(
            "row m/z" = 100,
            "sample1.mzML Feature status" = "DETECTED",
            check.names = FALSE
        ),
        feature_path,
        row.names = FALSE
    )

    # act, assert
    expect_error(
        MZmineFeatureListsSource(feature_path, sample_paths),
        "Missing required columns"
    )
})

test_that("MZmineFeatureListsSource reads an aligned MZmine >=3 feature list", {
    # arrange: two features aligned across two files; the second feature has
    # no peak in sampleB (UNKNOWN) and a gap-filled one in sampleA.
    tmp <- tempfile("mzmine")
    dir.create(tmp)
    sample_paths <- file.path(tmp, c("sampleA.mzML", "sampleB.mzML"))
    file.create(sample_paths)
    feature_path <- file.path(tmp, "aligned.csv")

    per_file <- function(file, mz, rt, state) {
        df <- data.frame(mz, rt, area = 1000, height = 100, mz - 0.01,
                         mz + 0.01, rt - 0.1, rt + 0.1, state)
        names(df) <- paste0("datafile:", file, ":", c(
            "mz", "rt", "area", "height", "mz_range:min", "mz_range:max",
            "rt_range:min", "rt_range:max", "feature_state"))
        df
    }
    write.csv(
        cbind(
            data.frame(id = 1:2, mz = c(150, 200), rt = c(10, 12),
                       check.names = FALSE),
            per_file("sampleA.mzML", c(150.001, 200.001), c(10.1, 12.1),
                     c("DETECTED", "ESTIMATED")),
            per_file("sampleB.mzML", c(149.999, NA), c(9.9, NA),
                     c("DETECTED", "UNKNOWN"))),
        feature_path, row.names = FALSE)

    # act
    res <- MZmineFeatureListsSource(feature_path, sample_paths)

    # assert: each file's own values, without the missing peak
    expect_equal(nrow(res@peaks), 3)
    expect_equal(sort(res@peaks$sample_index), c(1, 1, 2))
    a <- res@peaks[res@peaks$sample_index == 1, ]
    expect_equal(sort(a$mz), c(150.001, 200.001))
    expect_equal(sort(a$rt), c(10.1, 12.1) * 60)
    b <- res@peaks[res@peaks$sample_index == 2, ]
    expect_equal(b$mz, 149.999)
})

test_that("MZmineFeatureListsSource reads an aligned MZmine 2 feature list", {
    # arrange
    tmp <- tempfile("mzmine")
    dir.create(tmp)
    sample_paths <- file.path(tmp, c("s1.mzML", "s2.mzML"))
    file.create(sample_paths)
    feature_path <- file.path(tmp, "aligned.csv")

    per_file <- function(file, rt, state) {
        df <- data.frame(state, 100, rt, 1000, 500, 99.9, 100.1, rt - 0.2,
                         rt + 0.2)
        names(df) <- paste0(file, c(
            " Feature status", " Feature m/z", " Feature RT", " Peak area",
            " Peak height", " Feature m/z min", " Feature m/z max",
            " Feature RT start", " Feature RT end"))
        df
    }
    write.csv(
        cbind(data.frame("row m/z" = 100, check.names = FALSE),
              per_file("s1.mzML", 5, "DETECTED"),
              per_file("s2.mzML", 5.1, "DETECTED")),
        feature_path, row.names = FALSE)

    # act
    res <- MZmineFeatureListsSource(feature_path, sample_paths)

    # assert
    expect_equal(nrow(res@peaks), 2)
    expect_equal(res@peaks$rt[order(res@peaks$sample_index)], c(5, 5.1) * 60)
})
