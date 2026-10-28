export_pipeline <- function() {
  lcmsPlot(get_test_sample_paths()) +
    lp_chromatogram(features = rbind(c(
      mzmin = 334.9, mzmax = 335.1, rtmin = 2700, rtmax = 2900))) +
    lp_arrange(group_by = "sample_id")
}

new_temp_dir <- function() {
  dir <- tempfile("export-")
  dir.create(dir)
  dir
}

test_that("lp_save writes the figure and its provenance", {
  # arrange
  p <- export_pipeline()
  dir <- new_temp_dir()

  # act
  pdf_file <- lp_save(p, file.path(dir, "figure.pdf"))
  png_file <- lp_save(p, file.path(dir, "figure.png"), dpi = 72)

  # assert
  expect_true(file.exists(pdf_file))
  expect_true(file.exists(png_file))
  pdf_prov <- jsonlite::read_json(paste0(pdf_file, ".provenance.json"))
  png_prov <- jsonlite::read_json(paste0(png_file, ".provenance.json"))
  expect_equal(pdf_prov$output$format, "pdf")
  expect_equal(png_prov$output$format, "png")
  expect_equal(png_prov$output$dpi, 72)
  expect_null(pdf_prov$output$data)
  expect_false(dir.exists(file.path(dir, "figure_data")))
})

test_that("lp_save writes the plotted data with data = TRUE", {
  # arrange
  p <- export_pipeline()
  dir <- new_temp_dir()

  # act
  lp_save(p, file.path(dir, "figure.pdf"), data = TRUE)

  # assert
  expect_setequal(
    dir(file.path(dir, "figure_data")),
    c("chromatograms.csv", "samples.csv", "columns.csv"))
  prov <- jsonlite::read_json(file.path(dir, "figure.pdf.provenance.json"))
  expect_equal(prov$output$data, "figure_data")
})

test_that("lp_save skips the provenance with provenance = FALSE", {
  # arrange
  p <- export_pipeline()
  file <- tempfile(fileext = ".pdf")

  # act
  lp_save(p, file, provenance = FALSE)

  # assert
  expect_true(file.exists(file))
  expect_false(file.exists(paste0(file, ".provenance.json")))
})

test_that("lp_export_data writes the plotted datasets and a data dictionary", {
  # arrange
  p <- export_pipeline()
  dir <- new_temp_dir()

  # act
  paths <- lp_export_data(p, dir)

  # assert
  expect_equal(names(paths), c("chromatograms", "samples", "columns"))
  plotted <- .plot_datasets(p, .additional_datasets(p))$chromatograms
  written <- read.csv(paths[["chromatograms"]])
  expect_equal(nrow(written), nrow(plotted))
  expect_equal(written$rt, plotted$rt)
  expect_equal(written$intensity, plotted$intensity)
  expect_equal(read.csv(paths[["samples"]])$sample_id, sampleMetadata(p)$sample_id)

  columns <- read.csv(paths[["columns"]])
  expect_equal(
    columns$column[columns$table == "chromatograms"], colnames(plotted))
  rt <- columns[columns$table == "chromatograms" & columns$column == "rt", ]
  expect_equal(rt$unit, "s")
})

test_that("lp_export_data writes the detected peaks and tab-separated files", {
  # arrange
  p <- lcmsPlot(cached_xcms_object(), sample_id_column = "sample_name") +
    lp_chromatogram(features = baseline_feature_ids(), highlight_peaks = TRUE)
  dir <- new_temp_dir()

  # act
  paths <- lp_export_data(p, dir, format = "tsv")

  # assert
  expect_true(all(grepl("\\.tsv$", paths)))
  peaks <- read.delim(paths[["detected_peaks"]])
  expect_equal(nrow(peaks), nrow(detectedPeaks(p)))
})

test_that("lp_save and lp_export_data need data to plot", {
  # arrange
  p <- lcmsPlot(get_test_sample_paths())

  # act / assert
  expect_error(lp_save(p, tempfile(fileext = ".pdf")), "lp_save")
  expect_error(lp_export_data(p, tempdir()), "lp_export_data")
})
