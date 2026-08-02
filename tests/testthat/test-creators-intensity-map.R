test_that("create_intensity_map survives a window matching one peak per scan", {
  # arrange: a +/- 0.5 Da window around a real mass trace leaves most scans with
  # exactly one in-range peak, which used to collapse the matrix to a vector
  raw_files <- get_test_sample_paths(1)

  # act
  obj <- lcmsPlot(raw_files) +
    lp_intensity_map(mz_range = c(495.7, 496.7), rt_range = c(2500, 3500))

  # assert
  expect_gt(nrow(obj@data@intensity_maps), 0)
})

test_that("create_intensity_map emits the columns the slot validator pins", {
  # arrange
  raw_files <- get_test_sample_paths(1)

  # act
  obj <- lcmsPlot(raw_files) +
    lp_intensity_map(mz_range = c(300, 320), rt_range = c(2500, 2700))

  # assert
  expect_equal(
    colnames(obj@data@intensity_maps),
    c("rt", "mz", "intensity", "metadata_index", "feature_metadata_id"))
})

test_that("create_intensity_map skips binning for the point geom", {
  # arrange
  raw_files <- get_test_sample_paths(1)

  # act
  tiled <- lcmsPlot(raw_files) +
    lp_intensity_map(mz_range = c(300, 320), rt_range = c(2500, 2700))
  pointed <- lcmsPlot(raw_files) +
    lp_intensity_map(mz_range = c(300, 320), rt_range = c(2500, 2700),
                     geom = "point")

  # assert
  expect_equal(
    colnames(pointed@data@intensity_maps),
    colnames(tiled@data@intensity_maps))
  # binned m/z values sit on the 0.1 grid; raw centroids do not
  expect_true(all(tiled@data@intensity_maps$mz ==
                    round(tiled@data@intensity_maps$mz, 1)))
  expect_false(all(pointed@data@intensity_maps$mz ==
                     round(pointed@data@intensity_maps$mz, 1)))
})

test_that("create_intensity_map honours the bin widths", {
  # arrange
  raw_files <- get_test_sample_paths(1)

  # act
  coarse <- lcmsPlot(raw_files) +
    lp_intensity_map(mz_range = c(300, 320), rt_range = c(2500, 2700),
                     bin_rt = 10, bin_mz = 1)

  # assert
  expect_true(all(coarse@data@intensity_maps$mz %% 1 == 0))
  expect_true(all(coarse@data@intensity_maps$rt %% 10 == 0))
})

test_that("bin_coordinate matches round() for power-of-ten widths", {
  # arrange
  set.seed(42)
  x <- stats::runif(1000, 0, 1000)

  # act / assert
  expect_identical(lcmsPlot:::bin_coordinate(x, 0.1), round(x, 1))
  expect_identical(lcmsPlot:::bin_coordinate(x, 0.01), round(x, 2))
  expect_identical(lcmsPlot:::bin_coordinate(x, 1), round(x, 0))
  expect_true(all(lcmsPlot:::bin_coordinate(x, 30) %% 30 == 0))
})

test_that("plot_intensity_map dispatches on the geom option", {
  # arrange
  dataset <- list(intensity_maps = data.frame(
    rt = c(100, 100, 200), mz = c(300.1, 300.2, 300.1),
    intensity = c(10, 20, 30),
    metadata_index = 1, feature_metadata_id = 1))
  opts <- lcmsPlot:::default_options()

  geom_of <- function(geom) {
    o <- opts
    o$intensity_maps$geom <- geom
    class(plot_intensity_map(dataset, NULL, o, single = TRUE)$layers[[1]]$geom)[1]
  }

  # act / assert
  expect_equal(geom_of("tile"), "GeomTile")
  expect_equal(geom_of("point"), "GeomPoint")
  expect_equal(geom_of("density"), "GeomDensity2dFilled")
})

test_that("lp_intensity_map records the geom and rejects unknown ones", {
  # arrange
  dataset <- list(intensity_maps = data.frame(
    rt = c(100, 100, 200), mz = c(300.1, 300.2, 300.1),
    intensity = c(10, 20, 30),
    metadata_index = 1, feature_metadata_id = 1))

  build <- function(...) {
    obj <- lcmsPlot(get_test_sample_paths(1))
    opts <- lp_intensity_map(mz_range = c(300, 320), rt_range = c(2500, 2510),
                             ...)(obj)@options
    plot_intensity_map(dataset, NULL, opts, single = TRUE)$layers[[1]]$geom
  }

  # act / assert
  expect_true(inherits(build(geom = "density"), "GeomDensity2dFilled"))
  expect_true(inherits(build(geom = "point"), "GeomPoint"))
  expect_true(inherits(build(), "GeomTile"))
  expect_error(
    lp_intensity_map(mz_range = c(300, 320), rt_range = c(2500, 2510),
                     geom = "hexbin"),
    "should be one of")
})
