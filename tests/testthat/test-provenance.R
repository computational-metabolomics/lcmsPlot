provenance_pipeline <- function() {
  features <- data.frame(
    feature_id = "A", mzmin = 334.9, mzmax = 335.1, rtmin = 2700, rtmax = 2900)
  lcmsPlot(get_test_sample_paths()) +
    lp_chromatogram(features = features) +
    lp_arrange(group_by = "sample_id")
}

test_that(".provenance records the versions, the pipeline and the inputs", {
  # arrange
  p <- provenance_pipeline()

  # act
  prov <- .provenance(p, output = list(file = "figure.pdf"))

  # assert
  expect_equal(prov$generator, paste("lcmsPlot", packageVersion("lcmsPlot")))
  expect_match(prov$generated_at, "^\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}:\\d{2}Z$")
  expect_equal(prov$output$file, "figure.pdf")
  expect_equal(prov$environment$r, R.version$version.string)
  expect_true("ggplot2" %in% names(prov$environment$packages))
  expect_equal(
    vapply(prov$pipeline$layers, `[[`, character(1), "layer"),
    c("lp_chromatogram", "lp_arrange"))
  expect_equal(prov$pipeline$layers[[1]]$args$features$feature_id, "A")
  expect_equal(prov$input$source, "character")
  expect_equal(
    prov$input$files$path, sampleMetadata(p)$sample_path)
  expect_equal(
    prov$input$files$size, unname(file.size(get_test_sample_paths())))
  expect_false("md5" %in% colnames(prov$input$files))
})

test_that(".provenance records the MD5 checksums only on request", {
  # arrange
  p <- provenance_pipeline()

  # act
  prov <- .provenance(p, checksums = TRUE)

  # assert
  expect_equal(
    prov$input$files$md5, unname(tools::md5sum(get_test_sample_paths())))
})

test_that(".json_safe writes functions, matrices and objects", {
  # act
  safe <- .json_safe(list(
    fun = function(x) max(x),
    matrix = rbind(c(mzmin = 1, mzmax = 2)),
    factor = factor("a"),
    object = BiocParallel::SerialParam()))

  # assert
  expect_match(safe$fun, "max(x)", fixed = TRUE)
  expect_equal(safe$matrix, data.frame(mzmin = 1, mzmax = 2))
  expect_equal(safe$factor, "a")
  expect_equal(safe$object, list(class = "SerialParam"))
})

test_that("the provenance round-trips through JSON", {
  # arrange
  p <- provenance_pipeline()
  file <- tempfile(fileext = ".pdf")

  # act
  path <- .write_provenance(.provenance(p), file)

  # assert
  expect_equal(path, paste0(file, ".provenance.json"))
  prov <- jsonlite::read_json(path, simplifyVector = TRUE)
  expect_equal(prov$pipeline$layers$layer, c("lp_chromatogram", "lp_arrange"))
  expect_equal(prov$input$files$sample_id, sampleMetadata(p)$sample_id)
})
