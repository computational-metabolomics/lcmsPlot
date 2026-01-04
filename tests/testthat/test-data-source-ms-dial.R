test_that("MsDialPeaksSource creates ExternalDataSource", {
    # arrange
    tmp <- tempdir()

    sample_paths <- file.path(tmp, "sample1.mzML")
    file.create(sample_paths)

    peak_path <- file.path(tmp, "peaks.csv")

    write.csv(
        data.frame(
            "Precursor m/z" = 300,
            "RT (min)" = 6,
            "Area" = 1500,
            "Height" = 700,
            "RT left(min)" = 5.8,
            "RT right (min)" = 6.2,
            check.names = FALSE
        ),
        peak_path,
        row.names = FALSE
    )

    # act
    res <- MsDialPeaksSource(
        peaks_paths = peak_path,
        sample_paths = sample_paths
    )

    # assert
    expect_s4_class(res, "ExternalDataSource")
    expect_equal(res@name, "ms-dial")

    expect_true(is.data.frame(res@metadata))
    expect_true(is.data.frame(res@peaks))

    expect_equal(nrow(res@peaks), 1)

    expect_true(all(
        c("mz", "rt", "into", "maxo", "rtmin", "rtmax", "sample_index") %in%
            colnames(res@peaks)
    ))

    expect_equal(res@peaks$sample_index, 1)
})

test_that("Multiple MS-DIAL peak lists are combined correctly", {
    # arrange
    tmp <- tempdir()

    sample_paths <- file.path(tmp, c("s1.mzML", "s2.mzML"))
    file.create(sample_paths)

    make_peaks <- function(mz) {
        data.frame(
            "Precursor m/z" = mz,
            "RT (min)" = 5,
            "Area" = 1000,
            "Height" = 400,
            "RT left(min)" = 4.8,
            "RT right (min)" = 5.2,
            check.names = FALSE
        )
    }

    paths <- file.path(tmp, c("p1.csv", "p2.csv"))
    write.csv(make_peaks(100), paths[1], row.names = FALSE)
    write.csv(make_peaks(200), paths[2], row.names = FALSE)

    # act
    res <- MsDialPeaksSource(paths, sample_paths)

    # assert
    expect_equal(nrow(res@peaks), 2)
    expect_equal(res@peaks$sample_index, c(1, 2))
})

test_that("MsDialPeaksSource fails when required MS-DIAL columns are missing", {
    # arrange
    tmp <- tempdir()

    sample_paths <- file.path(tmp, "sample1.mzML")
    file.create(sample_paths)

    peak_path <- file.path(tmp, "peaks.csv")

    write.csv(
        data.frame(
            "Precursor m/z" = 300,
            "RT (min)" = 6,
            check.names = FALSE
        ),
        peak_path,
        row.names = FALSE
    )

    # act, assert
    expect_error(
        MsDialPeaksSource(peak_path, sample_paths),
        "Missing required columns in MS-DIAL peak list"
    )
})
