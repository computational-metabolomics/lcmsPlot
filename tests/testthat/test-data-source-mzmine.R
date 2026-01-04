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
