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

highlight_peaks <- function(data, options) {
  if (options$chromatograms$highlight_peaks) {
    return(geom_ribbon(
      data = subset(data, rt >= peak_rt_min & rt <= peak_rt_max),
      mapping = highlight_peaks_aes(ymax = "intensity", options),
      alpha = 0.3
    ))
  } else {
    return(NULL)
  }
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
    return(facet_wrap(as.formula(paste("~", paste(options$facets$facets, collapse = "+")))))
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
