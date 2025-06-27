build_aes <- function(x, y, options) {
  p_aes <- list(x = sym(x), y = sym(y))

  if (!is.null(options$arrangement$group_by)) {
    p_aes$color <- sym(options$arrangement$group_by)
  }

  return(do.call(aes, p_aes))
}

highlight_peaks_aes <- function(ymax, options) {
  p_aes <- list(ymax = sym(ymax), ymin = 0)

  if (!is.null(options$arrangement$group_by)) {
    p_aes$fill <- sym(options$arrangement$group_by)
  }

  return(do.call(aes, p_aes))
}

highlight_peaks <- function(dataset, detected_peaks, options) {
  if (options$chromatograms$highlight_peaks) {
    highlight_df <- detected_peaks %>%
      group_by(sample_id) %>%
      mutate(peak_id = row_number()) %>%
      ungroup() %>%
      rowwise() %>%
      do({
        peak <- .
        subset <- dataset %>%
          filter(sample_id == peak$sample_id) %>%
          filter(rt >= peak$rtmin, rt <= peak$rtmax) %>%
          mutate(peak_id = peak$peak_id)
        subset
      }) %>%
      bind_rows()

    return(geom_ribbon(
      data = highlight_df,
      aes(
        ymin = 0,
        ymax = intensity,
        group = peak_id,
        fill = factor(sample_id),
        colour = factor(sample_id)), # TODO: check this, should it be peak or sample?
      alpha = 0.3,
      linetype = 1
    ))
  } else {
    return(NULL)
  }
}

highlight_spectra_scans <- function(dataset, options) {
  geom_vline(data = dataset, aes(xintercept = rt), color = "black", linetype = "dashed")
}

rt_lines <- function(options) {
  lapply(options$rt_lines, function (rt_line_obj) {
    geom_vline(xintercept = rt_line_obj$intercept, color = rt_line_obj$color, linetype = rt_line_obj$line_type)
  })
}

legend_title <- function(options) {
  if (!is.null(options$labels$legend)) {
    return(guides(
      color = guide_legend(title = options$labels$legend),
      fill = guide_legend(title = options$labels$legend)
    ))
  } else {
    return(NULL)
  }
}

faceting <- function(options, single) {
  if (single & !is.null(options$facets$facets)) {
    return(facet_wrap(
      as.formula(paste("~", paste(options$facets$facets, collapse = "+"))),
      ncol = options$facets$ncol,
      nrow = options$facets$nrow))
  } else {
    return(NULL)
  }
}

grid_layout <- function(options, single) {
  rows <- options$grid$rows
  cols <- options$grid$cols

  if (single & (!is.null(rows) | !is.null(cols))) {
    facet_formula <- if (!is.null(rows) & !is.null(cols)) {
      as.formula(paste(rows, "~", cols))
    } else if (!is.null(rows)) {
      as.formula(paste(rows, "~ ."))
    } else if (!is.null(cols)) {
      as.formula(paste(". ~", cols))
    } else {
      NULL  # No faceting
    }

    return(facet_grid(facet_formula))
  } else {
    return(NULL)
  }
}

title_tile <- function(text) {
  text_grob <- grid::textGrob(text, gp = grid::gpar(col = "black", fontsize = 9))
  return(text_grob)
  # gray_rect <- grid::rectGrob(gp = grid::gpar(fill = "gray90", col = NA))
  # grid::grobTree(gray_rect, text_grob)
}
