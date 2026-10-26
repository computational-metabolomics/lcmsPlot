theme_opts <- function(group_by = "sample_id", ...) {
  opts <- lcmsPlot:::default_options()
  opts$arrangement$group_by <- group_by
  opts$theme <- utils::modifyList(opts$theme, list(...))
  opts
}

two_group_chromatograms <- function() {
  tibble::tibble(
    rt = rep(1:5, 4),
    intensity = rep(c(1, 4, 9, 4, 1), 4),
    sample_id = rep(c("s1", "s2", "s3", "s4"), each = 5),
    sample_group = rep(c("KO", "KO", "WT", "WT"), each = 5))
}

test_that(".palette_colours picks Okabe-Ito and falls back to viridis", {
  expect_equal(.palette_colours("okabe_ito", 3), .okabe_ito[1:3])
  expect_equal(.palette_colours("okabe_ito", 8), .okabe_ito)
  expect_equal(.palette_colours("okabe_ito", 9), scales::viridis_pal()(9))
  expect_equal(.palette_colours("viridis", 4), scales::viridis_pal()(4))
  expect_equal(.palette_colours("ggplot2", 4), scales::hue_pal()(4))
})

test_that(".palette_colours uses a custom palette in order", {
  expect_equal(.palette_colours(c("red", "blue", "green"), 2), c("red", "blue"))
  expect_error(.palette_colours(c("red", "blue"), 3), "has 2 colours")
})

test_that("lp_theme validates its arguments", {
  expect_error(lp_theme(preset = "fancy"), "'preset'")
  expect_error(lp_theme(palette = c("red", "not-a-colour")), "not-a-colour")
  expect_error(lp_theme(palette = 1:3), "'palette'")
  expect_error(lp_theme(base_size = -1), "'base_size'")
  expect_error(lp_theme(line_width = c(1, 2)), "'line_width'")
  expect_no_error(lp_theme(preset = ggplot2::theme_bw()))
  expect_no_error(lp_theme(palette = c(KO = "#0072B2", WT = "#D55E00")))
})

test_that("lp_theme stores its options", {
  obj <- lcmsPlot(get_test_sample_paths()) +
    lp_theme(preset = "publication", palette = "viridis", base_size = 8)

  expect_equal(obj@options$theme$preset, "publication")
  expect_equal(obj@options$theme$palette, "viridis")
  expect_equal(obj@options$theme$base_size, 8)
  expect_null(obj@options$theme$line_width)
  expect_equal(obj@history[[1]]$name, "lp_theme")
})

test_that(".base_theme builds each preset", {
  minimal <- .base_theme(theme_opts())
  expect_equal(minimal$text$size, 11)

  publication <- .base_theme(theme_opts(preset = "publication"))
  expect_equal(publication$text$size, 8)
  expect_s3_class(publication$panel.grid.minor, "element_blank")
  expect_false(inherits(publication$axis.line, "element_blank"))
  expect_equal(publication$panel.background$fill, "white")

  sized <- .base_theme(theme_opts(preset = "classic", base_size = 14))
  expect_equal(sized$text$size, 14)

  custom <- ggplot2::theme_dark()
  expect_identical(.base_theme(theme_opts(preset = custom)), custom)
})

test_that(".series_values assigns colours in ggplot2's order", {
  datasets <- list(chromatograms = two_group_chromatograms())

  values <- .series_values(datasets, NULL, theme_opts("sample_group"))
  expect_equal(values, c(KO = .okabe_ito[1], WT = .okabe_ito[2]))

  datasets$chromatograms$sample_group <- factor(
    datasets$chromatograms$sample_group, levels = c("WT", "KO"))
  values <- .series_values(datasets, NULL, theme_opts("sample_group"))
  expect_equal(names(values), c("WT", "KO"))
})

test_that(".series_values keeps a named palette and skips continuous columns", {
  datasets <- list(chromatograms = two_group_chromatograms())
  named <- c(KO = "black", WT = "grey")

  expect_equal(
    .series_values(datasets, NULL, theme_opts("sample_group", palette = named)),
    named)
  expect_null(.series_values(datasets, NULL, theme_opts("rt")))
  expect_null(.series_values(datasets, NULL, theme_opts(NULL)))
})

test_that("chromatogram traces take the palette and the line width", {
  datasets <- list(chromatograms = two_group_chromatograms())
  opts <- theme_opts(line_width = 1.2)
  opts$theme$series_values <- .series_values(datasets, NULL, opts)

  p <- plot_chromatogram(datasets, list(detected_peaks = tibble::tibble()), opts)
  data <- ggplot2::layer_data(p, 1)

  expect_setequal(unique(data$colour), .okabe_ito[1:4])
  expect_equal(unique(data$linewidth), 1.2)
})

test_that("more than eight series fall back to viridis", {
  chromatograms <- tibble::tibble(
    rt = rep(1:3, 9),
    intensity = 1,
    sample_id = rep(sprintf("s%02d", 1:9), each = 3))
  datasets <- list(chromatograms = chromatograms)
  opts <- theme_opts()
  opts$theme$series_values <- .series_values(datasets, NULL, opts)

  p <- plot_chromatogram(datasets, list(detected_peaks = tibble::tibble()), opts)

  expect_setequal(
    unique(ggplot2::layer_data(p, 1)$colour), scales::viridis_pal()(9))
})

test_that("highlighted peaks take the colour of their trace", {
  datasets <- list(chromatograms = two_group_chromatograms())
  peaks <- tibble::tibble(
    sample_id = c("s1", "s3"), rt = 3, rtmin = 2, rtmax = 4, maxo = 9)

  for (group_by in c("sample_id", "sample_group")) {
    opts <- theme_opts(group_by)
    opts$chromatograms$highlight_peaks <- TRUE
    opts$theme$series_values <- .series_values(datasets, peaks, opts)

    p <- plot_chromatogram(datasets, list(detected_peaks = peaks), opts)
    traces <- ggplot2::layer_data(p, 1)
    ribbons <- ggplot2::layer_data(p, 2)

    # One trace per sample, grouped in sample order
    trace_colours <- unique(traces[, c("group", "colour")])$colour
    ribbon_rows <- ribbons[!duplicated(ribbons$group), ]

    expect_equal(ribbon_rows$fill, trace_colours[c(1, 3)], info = group_by)
    expect_equal(ribbon_rows$colour, ribbon_rows$fill, info = group_by)
  }
})

test_that("an explicit highlight colour is left alone", {
  datasets <- list(chromatograms = two_group_chromatograms())
  peaks <- tibble::tibble(
    sample_id = "s1", rt = 3, rtmin = 2, rtmax = 4, maxo = 9)
  opts <- theme_opts()
  opts$chromatograms$highlight_peaks <- TRUE
  opts$chromatograms$highlight_peaks_color <- "#121212"
  opts$theme$series_values <- .series_values(datasets, peaks, opts)

  p <- plot_chromatogram(datasets, list(detected_peaks = peaks), opts)

  expect_equal(unique(ggplot2::layer_data(p, 2)$fill), "#121212")
})

test_that("continuous scales are unchanged", {
  maps <- tibble::tibble(
    rt = rep(1:3, 3), mz = rep(1:3, each = 3), intensity = 1:9)
  opts <- theme_opts(preset = "publication")

  p <- plot_intensity_map(list(intensity_maps = maps), list(), opts)

  expect_s3_class(p$scales$get_scales("fill"), "ScaleContinuous")
  expect_equal(
    ggplot2::layer_data(p, 1)$fill[1],
    ggplot2::layer_data(
      ggplot2::ggplot(maps, ggplot2::aes(rt, mz, fill = log1p(intensity))) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_viridis_c(), 1)$fill[1])
})

test_that("the purity overlay keeps its gradient next to the series colours", {
  skip_if_not_installed("ggnewscale")

  datasets <- list(chromatograms = two_group_chromatograms())
  opts <- theme_opts()
  opts$theme$series_values <- .series_values(datasets, NULL, opts)
  opts$purity_overlay$show <- TRUE
  scores <- tibble::tibble(rt = c(2, 4), in_purity = c(0, 1), sample_index = 1L)

  p <- plot_chromatogram(
    datasets,
    list(detected_peaks = tibble::tibble(), purity_scores = scores),
    opts)
  overlay <- ggplot2::layer_data(p, length(p$layers))

  expect_equal(toupper(overlay$colour), c("#D73027", "#1A9850"))
  expect_setequal(unique(ggplot2::layer_data(p, 1)$colour), .okabe_ito[1:4])
})

test_that("every panel of a faceted layout shares the theme and the colours", {
  datasets <- list(
    chromatograms = tibble::tibble(
      rt = rep(1:3, 3), intensity = 1,
      feature_id = c(rep("A", 6), rep("B", 3)),
      sample_id = c(rep(c("s1", "s2"), each = 3), rep("s2", 3))),
    spectra = tibble::tibble(
      mz = c(100, 200), intensity = 100, rt = 2, reference = FALSE,
      feature_id = c("A", "B"), sample_id = "s2"))
  obj <- lcmsPlot(get_test_sample_paths()) +
    lp_facets(facets = "feature_id") +
    lp_arrange(group_by = "sample_id") +
    lp_theme(preset = "publication")

  p <- plot_data(datasets, obj)
  panels <- list(p[[1]][[1]], p[[1]][[2]], p[[2]][[1]], p[[2]][[2]])

  for (panel in panels) {
    expect_s3_class(panel$theme$panel.grid.minor, "element_blank")
    expect_equal(panel$theme$text$size, 8)
  }
  # s2 is the second series of panel A but the only series of panel B
  s2_colour <- .okabe_ito[2]
  expect_true(s2_colour %in% ggplot2::layer_data(panels[[1]], 1)$colour)
  expect_equal(unique(ggplot2::layer_data(panels[[3]], 1)$colour), s2_colour)
  expect_equal(unique(ggplot2::layer_data(panels[[4]], 1)$colour), s2_colour)
})

test_that("a theme added after lp_get_plot() overrides the preset", {
  datasets <- list(chromatograms = two_group_chromatograms())
  obj <- lcmsPlot(get_test_sample_paths()) +
    lp_theme(preset = "publication")

  p <- plot_data(datasets, obj) +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "red"))

  expect_equal(p[[1]]$theme$panel.background$fill, "red")
  expect_s3_class(p[[1]]$theme$panel.grid.minor, "element_blank")
})
