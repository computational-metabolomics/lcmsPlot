# The per-peak filter that build_ribbon_data() replaced, kept as a reference.
naive_ribbon_data <- function(dataset, detected_peaks) {
  rows <- lapply(seq_len(nrow(detected_peaks)), function(i) {
    peak <- detected_peaks[i, ]
    out <- dataset[dataset$sample_id == peak$sample_id &
                     dataset$rt >= peak$rtmin &
                     dataset$rt <= peak$rtmax, ]
    out$peak_id <- rep(i, nrow(out))
    out
  })
  do.call(rbind, rows)
}

test_that("build_ribbon_data matches the per-peak definition", {
  # arrange
  dataset <- tibble::tibble(
    sample_id = rep(c("a", "b"), each = 11),
    rt = rep(0:10, 2),
    intensity_plot = seq_len(22),
    feature_id = rep(c("F1", "F2"), 11))
  detected_peaks <- tibble::tibble(
    sample_id = c("a", "a", "b", "c", "b"),
    # overlapping peaks in "a", boundary values, a sample with no rows ("c")
    rtmin = c(2, 4, 0, 1, 10),
    rtmax = c(5, 6, 3, 9, 10),
    mz = c(100, 101, 102, 103, 104))

  # act
  actual <- build_ribbon_data(dataset, detected_peaks)

  # assert
  expect_equal(
    as.data.frame(actual),
    as.data.frame(naive_ribbon_data(dataset, detected_peaks)))
  expect_named(actual, c(names(dataset), "peak_id"))
})

test_that("build_ribbon_data returns no rows when nothing falls under a peak", {
  # arrange
  dataset <- tibble::tibble(sample_id = "a", rt = 1:3, intensity_plot = 1:3)
  detected_peaks <- tibble::tibble(sample_id = "a", rtmin = 10, rtmax = 20)

  # act
  actual <- build_ribbon_data(dataset, detected_peaks)

  # assert
  expect_equal(nrow(actual), 0)
  expect_named(actual, c(names(dataset), "peak_id"))
})
