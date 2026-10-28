test_that("rt_diff holds the shift of every scan by the alignment", {
  # arrange
  p <- lcmsPlot(cached_xcms_object(), sample_id_column = "sample_name") +
    lp_rt_diff_plot()

  # act
  rt_diff <- summaries(p)$rt_diff

  # assert
  expect_gt(nrow(rt_diff), 0)
  expect_equal(rt_diff$diff, rt_diff$rt_adj - rt_diff$rt_raw)
  expect_true(any(rt_diff$diff != 0))
  expect_setequal(unique(rt_diff$sample_index), c(1L, 2L))
})
