# The pipeline of the report tests: three xcms features with their peaks.
report_pipeline <- function() {
  lcmsPlot(cached_xcms_object(), sample_id_column = "sample_name") +
    lp_chromatogram(features = baseline_feature_ids(), highlight_peaks = TRUE) +
    lp_arrange(group_by = "sample_group")
}

# The DDA standards with one feature per standard (see test-dataset-spectra.R).
report_dda_features <- function() {
  data.frame(
    sample_id = c("l-proline-MS2", "l-kynurenine-MS2"),
    feature_id = c("proline", "kynurenine"),
    mz = c(116.0705, 209.0918),
    rt = c(426, 365))
}

# The number of pages of a PDF file, from the /Count of its page tree.
pdf_page_count <- function(path) {
  bytes <- readBin(path, "raw", file.size(path))
  text <- rawToChar(bytes[bytes != as.raw(0)])
  counts <- regmatches(text, gregexpr("/Count [0-9]+", text, useBytes = TRUE))
  max(as.integer(sub("/Count ", "", counts[[1]])))
}

count_matches <- function(text, pattern) {
  lengths(regmatches(text, gregexpr(pattern, text, fixed = TRUE)))
}

test_that("lp_report writes an HTML page per feature with an index", {
  skip_if_not_installed("htmltools")
  # arrange
  p <- report_pipeline()
  file <- tempfile(fileext = ".html")

  # act
  index <- lp_report(p, file)

  # assert
  html <- paste(readLines(file, warn = FALSE), collapse = "\n")
  expect_equal(index$feature_id, baseline_feature_ids())
  expect_true(all(index$status == "ok"))
  expect_true(all(index$n_samples == 2))
  expect_equal(count_matches(html, "<section class=\"lp-report-page\""), 3)
  expect_equal(count_matches(html, "href=\"#page-"), 3)
  expect_equal(count_matches(html, "src=\"data:image/png;base64,"), 3)
  for (id in baseline_feature_ids()) {
    expect_match(html, paste0(">", id, "</a>"), fixed = TRUE)
  }
  expect_no_match(html, "<script[^>]+src=|<link[^>]+href=")
})

test_that("lp_report writes a multi-page PDF of the index and the features", {
  skip_if_not_installed("gridExtra")
  # arrange
  p <- report_pipeline()
  file <- tempfile(fileext = ".pdf")

  # act
  lp_report(p, file, page_size = "letter", orientation = "portrait")

  # assert
  # The index, a page per feature, and the provenance.
  expect_equal(pdf_page_count(file), 1 + 3 + 1)
})

test_that("lp_report records its provenance in the HTML and next to it", {
  skip_if_not_installed("htmltools")
  # arrange
  p <- report_pipeline()
  file <- tempfile(fileext = ".html")

  # act
  lp_report(p, file, checksums = TRUE)

  # assert
  html <- paste(readLines(file, warn = FALSE), collapse = "\n")
  generator <- paste("lcmsPlot", packageVersion("lcmsPlot"))
  expect_match(
    html, paste0("<meta name=\"generator\" content=\"", generator, "\""),
    fixed = TRUE)
  expect_equal(count_matches(html, "<section id=\"lp-provenance\">"), 1)
  expect_match(html, ">lp_chromatogram</td>", fixed = TRUE)
  embedded <- regmatches(html, regexpr(
    "(?s)id=\"lp-provenance-data\">.*?</script>", html, perl = TRUE))
  embedded <- sub("^id=\"lp-provenance-data\">", "", sub("</script>$", "", embedded))
  prov <- jsonlite::read_json(paste0(file, ".provenance.json"))
  expect_equal(jsonlite::fromJSON(embedded)$generated_at, prov$generated_at)
  expect_equal(prov$output$format, "html")
  expect_equal(prov$output$features, 3)
  expect_equal(prov$output$pages, 3)
  expect_length(prov$input$files, 2)
  expect_false(is.null(prov$input$files[[1]]$md5))
})

test_that("lp_report leaves out the provenance with provenance = FALSE", {
  skip_if_not_installed("htmltools")
  # arrange
  p <- report_pipeline()
  file <- tempfile(fileext = ".html")

  # act
  lp_report(p, file, provenance = FALSE)

  # assert
  html <- paste(readLines(file, warn = FALSE), collapse = "\n")
  expect_no_match(html, "lp-provenance", fixed = TRUE)
  expect_no_match(html, "name=\"generator\"", fixed = TRUE)
  expect_false(file.exists(paste0(file, ".provenance.json")))
})

test_that("lp_report writes the provenance of a PDF report next to it", {
  skip_if_not_installed("gridExtra")
  # arrange
  p <- report_pipeline()
  file <- tempfile(fileext = ".pdf")

  # act
  lp_report(p, file, provenance = FALSE)
  pages_without <- pdf_page_count(file)
  lp_report(p, file, title = "My report")

  # assert
  expect_equal(pdf_page_count(file), pages_without + 1)
  expect_true(file.exists(paste0(file, ".provenance.json")))
  bytes <- readBin(file, "raw", file.size(file))
  expect_true(grepl(
    "/Title (My report)", rawToChar(bytes[bytes != as.raw(0)]),
    fixed = TRUE, useBytes = TRUE))
})

test_that("lp_report draws several features per page", {
  skip_if_not_installed("htmltools")
  # arrange
  p <- report_pipeline()
  file <- tempfile(fileext = ".html")

  # act
  index <- lp_report(p, file, features_per_page = 2)

  # assert
  html <- paste(readLines(file, warn = FALSE), collapse = "\n")
  expect_equal(index$page, c(1L, 1L, 2L))
  expect_equal(count_matches(html, "<section class=\"lp-report-page\""), 2)
  expect_equal(count_matches(html, "href=\"#page-1\""), 2)
  expect_match(html, paste0(
    "<h2>", baseline_feature_ids()[1], ", ", baseline_feature_ids()[2],
    "</h2>"), fixed = TRUE)
})

test_that("lp_report fits several features per page in the PDF", {
  skip_if_not_installed("gridExtra")
  # arrange
  p <- report_pipeline()
  file <- tempfile(fileext = ".pdf")

  # act
  lp_report(p, file, features_per_page = 2, ncol = 1, theme = "journal")

  # assert
  expect_equal(pdf_page_count(file), 1 + 2 + 1)
})

test_that("lp_report leaves out the tables whose section is off", {
  skip_if_not_installed("htmltools")
  # arrange
  p <- report_pipeline()
  annotations_only <- tempfile(fileext = ".html")
  plots_only <- tempfile(fileext = ".html")

  # act
  lp_report(p, annotations_only, sections = "annotations")
  lp_report(p, plots_only, sections = character(), provenance = FALSE)

  # assert
  html <- paste(readLines(annotations_only, warn = FALSE), collapse = "\n")
  expect_equal(count_matches(html, "<h3>Annotations</h3>"), 3)
  expect_no_match(html, "Detected peaks", fixed = TRUE)
  html <- paste(readLines(plots_only, warn = FALSE), collapse = "\n")
  expect_no_match(html, "<h3>", fixed = TRUE)
  expect_equal(count_matches(html, "src=\"data:image/png;base64,"), 3)
})

test_that("lp_report styles the report with its theme", {
  skip_if_not_installed("htmltools")
  # arrange
  p <- report_pipeline()
  plain <- tempfile(fileext = ".html")
  clean <- tempfile(fileext = ".html")

  # act
  lp_report(p, plain, features = baseline_feature_ids()[1])
  lp_report(p, clean, features = baseline_feature_ids()[1], theme = "clean")

  # assert
  plain <- paste(readLines(plain, warn = FALSE), collapse = "\n")
  clean <- paste(readLines(clean, warn = FALSE), collapse = "\n")
  expect_no_match(plain, "--lp-accent", fixed = TRUE)
  expect_match(clean, "--lp-accent", fixed = TRUE)
  expect_match(
    clean, "<span class=\"lp-status lp-status-ok\">ok</span>", fixed = TRUE)
  expect_error(
    lp_report(p, tempfile(fileext = ".html"), theme = "neon"), "'arg'")
})

test_that("lp_report lists the features without data as empty", {
  skip_if_not_installed("htmltools")
  # arrange
  p <- report_pipeline()
  file <- tempfile(fileext = ".html")
  features <- c(baseline_feature_ids()[2], "missing")

  # act
  index <- lp_report(p, file, features = features)

  # assert
  html <- paste(readLines(file, warn = FALSE), collapse = "\n")
  expect_equal(index$feature_id, features)
  expect_equal(index$status, c("ok", "no data"))
  expect_equal(count_matches(html, "<section class=\"lp-report-page\""), 1)
  expect_match(html, "<td>missing</td>", fixed = TRUE)
})

test_that("lp_report marks a feature with no data points as empty", {
  skip_if_not_installed("htmltools")
  # arrange
  p <- lcmsPlot(get_dda_sample_paths()) +
    lp_chromatogram(
      features = rbind(nothing = c(mz = 900, rt = 10)), rt_tol = 1)
  file <- tempfile(fileext = ".html")

  # act
  index <- lp_report(p, file)

  # assert
  expect_equal(index$status, "no data")
  expect_equal(index$mz, 900)
})

test_that("lp_report pages by name with the annotations of the compounds", {
  skip_if_not_installed("htmltools")
  # arrange
  p <- lcmsPlot(baseline_lipid_search_source()) +
    lp_lipid_search() +
    lp_chromatogram(highlight_peaks = TRUE)
  file <- tempfile(fileext = ".html")

  # act
  index <- lp_report(p, file, by = "name")

  # assert
  html <- paste(readLines(file, warn = FALSE), collapse = "\n")
  expect_equal(index$name, c("PC(16:0_18:1)+H", "AcCa(14:0)+H"))
  expect_equal(index$formula, c("C42 H83 O8 N1 P1", "C21 H42 O4 N1"))
  expect_true(all(c("adduct", "mz", "rt") %in% colnames(index)))
  expect_equal(count_matches(html, "<section class=\"lp-report-page\""), 2)
  expect_match(html, "C42 H83 O8 N1 P1", fixed = TRUE)
})

test_that("lp_report shows the retention times in the unit of the plots", {
  skip_if_not_installed("htmltools")
  # arrange
  p <- lcmsPlot(cached_xcms_object(), sample_id_column = "sample_name") +
    lp_chromatogram(
      features = baseline_feature_ids()[1],
      highlight_peaks = TRUE,
      rt_unit = "minute")
  file <- tempfile(fileext = ".html")

  # act
  index <- lp_report(p, file)

  # assert
  html <- paste(readLines(file, warn = FALSE), collapse = "\n")
  expect_match(html, "<th>rt (min)</th>", fixed = TRUE)
  expect_match(html, sprintf("<td>%.2f</td>", index$rt / 60), fixed = TRUE)
})

test_that("lp_report does not extract the data again", {
  skip_if_not_installed("htmltools")
  # arrange
  p <- report_pipeline()
  file <- tempfile(fileext = ".html")
  local_mocked_bindings(
    create_chromatograms = function(...) stop("extracted again"),
    create_spectra = function(...) stop("extracted again"))

  # act
  index <- lp_report(p, file)

  # assert
  expect_true(all(index$status == "ok"))
})

test_that("lp_report extracts the MS2 spectra once for the whole report", {
  skip_if_not_installed("htmltools")
  # arrange
  p <- lcmsPlot(get_dda_sample_paths()) +
    lp_chromatogram(features = report_dda_features(), rt_tol = 15)
  file <- tempfile(fileext = ".html")
  calls <- 0
  create <- create_spectra
  local_mocked_bindings(create_spectra = function(...) {
    calls <<- calls + 1
    create(...)
  })

  # act
  suppressWarnings(lp_report(p, file, sections = "ms2"))

  # assert
  html <- paste(readLines(file, warn = FALSE), collapse = "\n")
  expect_equal(calls, 1)
  expect_equal(count_matches(html, "src=\"data:image/png;base64,"), 2)
  expect_no_match(html, "No MS2 spectra match", fixed = TRUE)
})

test_that("lp_report skips the MS2 section when the pipeline plots it", {
  skip_if_not_installed("htmltools")
  # arrange
  p <- suppressWarnings(
    lcmsPlot(get_dda_sample_paths()) +
      lp_chromatogram(features = report_dda_features(), rt_tol = 15) +
      lp_spectra(mode = "precursor"))
  file <- tempfile(fileext = ".html")

  # act / assert
  expect_message(
    suppressWarnings(lp_report(p, file, sections = "ms2")),
    "already plots the MS2 spectra")
})

test_that("lp_report keeps going when a page fails to render", {
  skip_if_not_installed("htmltools")
  # arrange
  p <- report_pipeline()
  file <- tempfile(fileext = ".html")
  local_mocked_bindings(.render_plot = function(...) stop("boom"))

  # act
  expect_warning(index <- lp_report(p, file), "3 page\\(s\\) failed")

  # assert
  expect_true(all(index$status == "error"))
  expect_true(file.exists(file))
})

test_that("lp_report rejects unknown formats, columns, and featureless data", {
  # arrange
  p <- report_pipeline()
  tic <- lcmsPlot(get_test_sample_paths()) +
    lp_chromatogram(aggregation_fun = "sum")

  # act / assert
  expect_error(lp_report(p, "report.docx"), "unknown format 'docx'")
  expect_error(
    lp_report(p, tempfile(fileext = ".html"), by = "nope"),
    "'by' must be one of")
  expect_error(
    lp_report(tic, tempfile(fileext = ".html")), "the data has no features")
  expect_error(
    lp_report(p, tempfile(fileext = ".html"), features_per_page = 0),
    "'features_per_page' must be")
  expect_error(
    lp_report(lcmsPlot(get_test_sample_paths()), "report.html"),
    "must be an lcmsPlot object with data")
})
