#' Render a built cell plot with base graphics
#'
#' Draws a projected 2D cell plot using the same resolved colors, sizes, alpha,
#' panel grid, theme, and annotations as the ggplot renderer. Animation frames
#' may supply shared axis limits.
#'
#' @param object A `cell_plot_built` object.
#' @param limits Optional shared x and y ranges used by animation frames.
#'
#' @return `NULL`, invisibly.
#'
#' @noRd
.render_cell_plot_base <- function(object, limits = NULL) {
  assert_class(object, "cell_plot_built", arg = "object")

  mapping <- object$mapping
  plot_data <- object$data
  n_nodes <- nrow(plot_data)
  colors <- scales::alpha(
    rep_len(.cell_plot_rendered_colors(object), n_nodes),
    alpha = rep_len(object$alpha$resolved, n_nodes)
  )
  cex <- .cell_relative_size_to_cex(
    .cell_plot_projected_sizes(object)
  )
  if (is.null(limits)) {
    limits <- .cell_animation_plot_limits(
      list(
        x = range(plot_data[[mapping$x]], na.rm = TRUE),
        y = range(plot_data[[mapping$y]], na.rm = TRUE)
      )
    )
  }

  row_levels <- if (is.null(object$grid$rows)) {
    list(NULL)
  } else {
    .cell_plot_facet_levels(plot_data[[object$grid$rows]])
  }
  col_levels <- if (is.null(object$grid$cols)) {
    list(NULL)
  } else {
    .cell_plot_facet_levels(plot_data[[object$grid$cols]])
  }
  n_row <- length(row_levels)
  n_col <- length(col_levels)
  has_legend <- !is.null(mapping$color)
  has_title <- !is.null(object$annotation$title)
  has_subtitle <- !is.null(object$annotation$subtitle)
  has_row_strips <- !is.null(object$grid$rows)
  has_col_strips <- !is.null(object$grid$cols)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(
    bg = object$theme$background_color,
    fg = object$theme$text_color,
    col.main = object$theme$text_color,
    col.lab = object$theme$text_color
  )

  text_cex <- object$theme$text_size / 12
  legend <- NULL
  if (has_legend) {
    legend <- .cell_base_legend_metrics(
      object,
      text_cex = text_cex,
      device_width = graphics::par("din")[[1]]
    )
  }

  layout <- .cell_base_layout_matrix(
    n_row = n_row,
    n_col = n_col,
    has_title = has_title,
    has_subtitle = has_subtitle,
    has_row_strips = has_row_strips,
    has_col_strips = has_col_strips,
    has_legend = has_legend,
    legend_width = legend$width
  )
  graphics::layout(
    layout$mat,
    widths = layout$widths,
    heights = layout$heights,
    respect = layout$respect
  )
  # layout() shrinks the base cex for grids with two or more rows or columns.
  # Reset it so text sizes follow the theme instead of the panel count.
  graphics::par(cex = 1)
  strip_cex <- text_cex * .cell_base_strip_text_scale

  if (has_title) {
    .cell_base_label(
      object$annotation$title,
      cex = text_cex * 1.2,
      color = object$theme$text_color,
      background = object$theme$background_color
    )
  }
  if (has_subtitle) {
    .cell_base_label(
      object$annotation$subtitle,
      cex = text_cex,
      color = object$theme$text_color,
      background = object$theme$background_color
    )
  }
  if (has_col_strips) {
    if (has_row_strips) {
      .cell_base_empty_panel(object$theme$background_color)
    }
    for (col_level in col_levels) {
      .cell_base_label(
        .cell_plot_facet_label(col_level),
        cex = strip_cex,
        color = object$theme$text_color,
        background = object$theme$strip_background_color
      )
    }
    if (has_legend) {
      .cell_base_empty_panel(object$theme$background_color)
    }
  }

  for (row_index in seq_len(n_row)) {
    if (has_row_strips) {
      .cell_base_label(
        .cell_plot_facet_label(row_levels[[row_index]]),
        cex = strip_cex,
        color = object$theme$text_color,
        background = object$theme$strip_background_color,
        srt = .cell_plot_row_strip_angle
      )
    }
    for (col_index in seq_len(n_col)) {
      keep <- rep(TRUE, n_nodes)
      if (has_row_strips) {
        keep <- keep &
          .cell_plot_facet_match(
            plot_data[[object$grid$rows]],
            row_levels[[row_index]]
          )
      }
      if (has_col_strips) {
        keep <- keep &
          .cell_plot_facet_match(
            plot_data[[object$grid$cols]],
            col_levels[[col_index]]
          )
      }
      .cell_base_scatter(
        x = plot_data[[mapping$x]][keep],
        y = plot_data[[mapping$y]][keep],
        colors = colors[keep],
        cex = if (length(cex) == 1L) cex else cex[keep],
        limits = limits,
        background = object$theme$background_color
      )
    }
  }

  if (has_legend) {
    .cell_base_legend(object, legend)
  }

  return(invisible(NULL))
}

#' Build a base-graphics layout for a cell plot
#'
#' @param n_row,n_col Panel counts.
#' @param has_title,has_subtitle,has_row_strips,has_col_strips,has_legend
#' Layout flags.
#' @param legend_width Legend column width in inches. Required when
#' `has_legend` is `TRUE`. The legend column is absolute so that the legend
#' title, color bar, and labels keep their physical size while the panels
#' share the remaining device width.
#'
#' @return A list with `mat`, `widths`, `heights`, and `respect`. `respect` is
#' `TRUE` only when both row and column strips are present, so that the
#' left-strip width equals the top-strip height in device units.
#'
#' @noRd
.cell_base_layout_matrix <- function(
  n_row,
  n_col,
  has_title,
  has_subtitle,
  has_row_strips,
  has_col_strips,
  has_legend,
  legend_width = NULL
) {
  if (has_legend) {
    assert_single_value(legend_width, type = "numeric", arg = "legend_width")
  }
  n_layout_row <- as.integer(has_title) + as.integer(has_subtitle) +
    as.integer(has_col_strips) + n_row
  n_layout_col <- as.integer(has_row_strips) + n_col + as.integer(has_legend)
  mat <- matrix(0L, nrow = n_layout_row, ncol = n_layout_col)
  id <- 1L
  layout_row <- 1L

  if (has_title) {
    mat[layout_row, ] <- id
    id <- id + 1L
    layout_row <- layout_row + 1L
  }
  if (has_subtitle) {
    mat[layout_row, ] <- id
    id <- id + 1L
    layout_row <- layout_row + 1L
  }
  if (has_col_strips) {
    start_col <- as.integer(has_row_strips) + 1L
    if (has_row_strips) {
      mat[layout_row, 1L] <- id
      id <- id + 1L
    }
    mat[layout_row, start_col:(start_col + n_col - 1L)] <-
      id + seq_len(n_col) - 1L
    id <- id + n_col
    if (has_legend) {
      mat[layout_row, n_layout_col] <- id
      id <- id + 1L
    }
    layout_row <- layout_row + 1L
  }

  panel_start_row <- layout_row
  for (row_index in seq_len(n_row)) {
    col_index <- 1L
    if (has_row_strips) {
      mat[layout_row, col_index] <- id
      id <- id + 1L
      col_index <- col_index + 1L
    }
    mat[layout_row, col_index:(col_index + n_col - 1L)] <-
      id + seq_len(n_col) - 1L
    id <- id + n_col
    layout_row <- layout_row + 1L
  }
  if (has_legend) {
    mat[panel_start_row:n_layout_row, n_layout_col] <- id
  }

  widths <- c(
    if (has_row_strips) 0.14,
    rep(1, n_col),
    if (has_legend) graphics::lcm(legend_width * 2.54)
  )
  heights <- c(
    if (has_title) 0.16,
    if (has_subtitle) 0.14,
    if (has_col_strips) 0.14,
    rep(1, n_row)
  )
  # Equal width and height units are only needed to make the left-strip
  # width match the top-strip height. Otherwise, let the layout fill the
  # device so titles and legends keep their relative size.
  return(list(
    mat = mat,
    widths = widths,
    heights = heights,
    respect = has_row_strips && has_col_strips
  ))
}

#' Draw one projected scatter panel
#'
#' @param x,y Point coordinates.
#' @param colors Point colors, including alpha.
#' @param cex Point sizes.
#' @param limits Shared axis ranges.
#' @param background Panel background color.
#'
#' @return `NULL`, invisibly.
#'
#' @noRd
.cell_base_scatter <- function(x, y, colors, cex, limits, background) {
  graphics::par(mar = c(0, 0, 0, 0), pty = "m", xaxs = "i", yaxs = "i")
  graphics::plot(
    0,
    0,
    xlim = limits$x,
    ylim = limits$y,
    type = "n",
    axes = FALSE,
    xlab = "",
    ylab = "",
    asp = 1
  )
  graphics::rect(
    limits$x[1],
    limits$y[1],
    limits$x[2],
    limits$y[2],
    col = background,
    border = NA
  )
  if (length(x) > 0L) {
    graphics::points(x, y, pch = 16, col = colors, cex = cex)
  }
  return(invisible(NULL))
}

#' Draw a text panel used for titles and facet strips
#'
#' @param label Text to draw.
#' @param cex Text size.
#' @param color Text color.
#' @param background Panel background color.
#' @param srt Text rotation in degrees.
#'
#' @return `NULL`, invisibly.
#'
#' @noRd
.cell_base_label <- function(label, cex, color, background, srt = 0) {
  graphics::par(
    mar = c(0, 0, 0, 0),
    xaxs = "i",
    yaxs = "i"
  )
  graphics::plot.new()
  graphics::rect(0, 0, 1, 1, col = background, border = NA)
  graphics::text(
    0.5,
    0.5,
    labels = label,
    cex = cex,
    col = color,
    srt = srt,
    xpd = TRUE
  )
  return(invisible(NULL))
}

#' Draw an empty layout cell
#'
#' @param background Background color.
#'
#' @return `NULL`, invisibly.
#'
#' @noRd
.cell_base_empty_panel <- function(background) {
  graphics::par(
    mar = c(0, 0, 0, 0),
    xaxs = "i",
    yaxs = "i"
  )
  graphics::plot.new()
  graphics::rect(0, 0, 1, 1, col = background, border = NA)
  return(invisible(NULL))
}

# Facet strip text is drawn at 80% of the theme text size, matching the
# ggplot2 `strip.text` default.
.cell_base_strip_text_scale <- 0.8

# Legend geometry expressed in multiples of the theme text size, following the
# ggplot2 defaults: 5.5 pt margins and key spacing, 1.2-line keys, legend text
# at 80% of the base size, and tick marks at 20% of the key size. The color bar
# is wider than the ggplot2 key so the gradient stays readable in small frames.
.cell_base_legend_geometry <- list(
  margin = 0.5,
  spacing = 0.5,
  key = 1.2,
  bar_width = 1.5,
  bar_height = 6,
  line_height = 1.25,
  text_scale = 0.8,
  tick_length = 0.2,
  max_width_fraction = 0.4
)

#' Measure the legend before the base-graphics layout is created
#'
#' Computes the legend column width in inches from the rendered text so the
#' legend title and labels always fit. The column grows with its contents up
#' to `max_width_fraction` of the device width, and it is always kept below
#' the device width so an absolute layout column cannot make `plot.new()`
#' fail. Titles and labels wrap onto several lines inside that limit. All
#' sizes are derived from `theme$text_size`, so the legend keeps the same
#' physical proportions as the ggplot2 renderer regardless of frame size.
#'
#' A graphics device must be open because text is measured with
#' [graphics::strwidth()].
#'
#' @param object A `cell_plot_built` object with a color mapping.
#' @param text_cex Theme text size as a `cex` multiplier.
#' @param device_width Device width in inches.
#'
#' @return A list describing the legend: `type`, `width` (inches), `unit`
#' (inches per text-size unit), `title_lines`, `title_cex`, `label_cex`,
#' `labels`, `label_lines` (wrapped label text, one character vector per
#' label), and, for continuous scales, `fills`, `breaks`, and `limits`, or,
#' for categorical scales, `colors`.
#'
#' @noRd
.cell_base_legend_metrics <- function(object, text_cex, device_width) {
  geometry <- .cell_base_legend_geometry
  unit <- object$theme$text_size / 72
  title_cex <- text_cex
  label_cex <- text_cex * geometry$text_scale
  margin <- geometry$margin * unit

  metrics <- list(
    type = object$color$type,
    unit = unit,
    title_cex = title_cex,
    label_cex = label_cex
  )
  if (identical(object$color$type, "categorical")) {
    metrics$labels <- as.character(object$color$limits)
    metrics$colors <- unname(object$color$colors[object$color$limits])
    key_width <- geometry$key * unit
  } else {
    legend_scale <- .cell_base_continuous_legend_scale(object$color)
    metrics$fills <- legend_scale$fills
    metrics$breaks <- legend_scale$breaks
    metrics$labels <- legend_scale$labels
    metrics$limits <- object$color$limits
    key_width <- geometry$bar_width * unit
  }

  # The color bar or keys, their margins, and the gap before the labels.
  # Labels wrap inside whatever remains of the column cap.
  chrome <- margin + key_width + geometry$spacing * unit + margin
  legend_limit <- .cell_base_legend_width_limit(
    device_width = device_width,
    preferred = geometry$max_width_fraction * device_width,
    chrome = chrome
  )
  label_max <- max(legend_limit - chrome, 0)
  metrics$label_lines <- lapply(metrics$labels, function(label) {
    .cell_base_wrap_text(label, max_width = label_max, cex = label_cex)
  })
  label_width <- .cell_base_widest_line(metrics$label_lines, cex = label_cex)
  body_width <- chrome + label_width

  title <- .cell_plot_legend_title(object)
  metrics$title_lines <- .cell_base_wrap_text(
    title,
    max_width = max(legend_limit - 2 * margin, 0),
    cex = title_cex
  )
  title_width <- .cell_base_widest_line(list(metrics$title_lines), cex = title_cex)

  metrics$width <- min(max(body_width, title_width + 2 * margin), legend_limit)
  return(metrics)
}

#' Cap a legend column so it cannot consume the device
#'
#' Prefers `preferred` (a fraction of the device) and widens that up to the
#' color-bar chrome when the bar itself needs more room. The result is always
#' strictly below the device width, because an absolute `layout()` column
#' wider than the device makes `plot.new()` fail.
#'
#' @param device_width Device width in inches.
#' @param preferred Preferred legend width in inches.
#' @param chrome Width of the legend bar or keys and their margins, in inches.
#'
#' @return Legend column limit in inches.
#'
#' @noRd
.cell_base_legend_width_limit <- function(device_width, preferred, chrome) {
  if (!is.finite(device_width) || device_width <= 0) {
    return(0)
  }
  hard_cap <- 0.9 * device_width
  limit <- min(max(preferred, chrome), hard_cap)
  return(max(limit, 0))
}

#' Widest rendered line among wrapped legend text
#'
#' @param line_groups A list of character vectors, one vector per label.
#' @param cex Text size as a `cex` multiplier.
#'
#' @return Width in inches, or `0` when there is no text.
#'
#' @noRd
.cell_base_widest_line <- function(line_groups, cex) {
  widths <- vapply(line_groups, function(lines) {
    if (length(lines) == 0L) {
      return(0)
    }
    max(.cell_base_text_width(lines, cex = cex))
  }, numeric(1))
  if (length(widths) == 0L) {
    return(0)
  }
  return(max(widths))
}

#' Measure rendered text width in inches
#'
#' @param labels Character labels.
#' @param cex Text size as a `cex` multiplier.
#'
#' @return Text widths in inches, one per label.
#'
#' @noRd
.cell_base_text_width <- function(labels, cex) {
  return(graphics::strwidth(labels, units = "inches", cex = cex))
}

#' Wrap text to a maximum rendered width
#'
#' Breaks text at whitespace first. A word that is wider than `max_width` on
#' its own is split between characters and each piece is placed on its own
#' line so every line fits.
#'
#' @param text A single character string.
#' @param max_width Maximum line width in inches.
#' @param cex Text size as a `cex` multiplier.
#'
#' @return A character vector with one element per line. Empty text yields an
#' empty vector.
#'
#' @noRd
.cell_base_wrap_text <- function(text, max_width, cex) {
  text <- trimws(as.character(text %||% ""))
  if (!nzchar(text)) {
    return(character())
  }
  words <- strsplit(text, "\\s+")[[1]]

  lines <- character()
  current <- ""
  for (word in words) {
    if (.cell_base_text_width(word, cex = cex) > max_width) {
      if (nzchar(current)) {
        lines <- c(lines, current)
        current <- ""
      }
      lines <- c(
        lines,
        .cell_base_break_word(word, max_width = max_width, cex = cex)
      )
      next
    }
    candidate <- if (nzchar(current)) paste(current, word) else word
    if (nzchar(current) && .cell_base_text_width(candidate, cex = cex) > max_width) {
      lines <- c(lines, current)
      current <- word
    } else {
      current <- candidate
    }
  }
  if (nzchar(current)) {
    lines <- c(lines, current)
  }
  return(lines)
}

#' Split one word so each piece fits a maximum rendered width
#'
#' @param word A single word without whitespace.
#' @param max_width Maximum piece width in inches.
#' @param cex Text size as a `cex` multiplier.
#'
#' @return A character vector of word pieces.
#'
#' @noRd
.cell_base_break_word <- function(word, max_width, cex) {
  glyphs <- strsplit(word, "")[[1]]
  pieces <- character()
  current <- ""
  for (glyph in glyphs) {
    candidate <- paste0(current, glyph)
    fits <- .cell_base_text_width(candidate, cex = cex) <= max_width
    if (nzchar(current) && !fits) {
      pieces <- c(pieces, current)
      current <- glyph
    } else {
      current <- candidate
    }
  }
  return(c(pieces, current))
}

#' Draw a categorical or continuous legend
#'
#' Draws the legend in a coordinate system measured in inches so the title,
#' keys, color bar, and labels match the sizes measured by
#' `.cell_base_legend_metrics()`. The legend block is left-aligned and
#' vertically centered in the legend column, like a ggplot2 legend placed to
#' the right of the panels.
#'
#' @param object A `cell_plot_built` object.
#' @param metrics Legend measurements from `.cell_base_legend_metrics()`.
#'
#' @return `NULL`, invisibly.
#'
#' @noRd
.cell_base_legend <- function(object, metrics) {
  geometry <- .cell_base_legend_geometry
  unit <- metrics$unit
  margin <- geometry$margin * unit
  spacing <- geometry$spacing * unit
  line_height <- geometry$line_height * unit
  text_color <- object$theme$text_color

  graphics::par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  graphics::plot.new()
  panel <- graphics::par("pin")
  graphics::plot.window(xlim = c(0, panel[[1]]), ylim = c(0, panel[[2]]))
  graphics::rect(
    0,
    0,
    panel[[1]],
    panel[[2]],
    col = object$theme$background_color,
    border = NA
  )

  title_height <- length(metrics$title_lines) * line_height
  title_gap <- if (title_height > 0) spacing else 0
  available <- panel[[2]] - 2 * margin - title_height - title_gap

  label_line <- geometry$line_height * geometry$text_scale * unit
  max_label_lines <- 1L
  if (length(metrics$label_lines) > 0L) {
    max_label_lines <- max(1L, lengths(metrics$label_lines))
  }
  if (identical(metrics$type, "categorical")) {
    key <- geometry$key * unit
    n_keys <- length(metrics$labels)
    row_height <- max(key, max_label_lines * label_line) + spacing
    if (n_keys > 0 && n_keys * row_height > available) {
      row_height <- max(available / n_keys, metrics$label_cex * 12 / 72)
    }
    body_height <- n_keys * row_height
  } else {
    body_height <- min(geometry$bar_height * unit, available)
    body_height <- max(body_height, 2 * unit)
  }

  block_height <- title_height + title_gap + body_height
  top <- min((panel[[2]] + block_height) / 2, panel[[2]] - margin)
  left <- margin

  if (title_height > 0) {
    graphics::text(
      left,
      top - (seq_along(metrics$title_lines) - 0.5) * line_height,
      labels = metrics$title_lines,
      adj = c(0, 0.5),
      cex = metrics$title_cex,
      col = text_color,
      xpd = TRUE
    )
  }
  body_top <- top - title_height - title_gap

  if (identical(metrics$type, "categorical")) {
    centers <- body_top - (seq_len(n_keys) - 0.5) * row_height
    graphics::points(
      rep(left + key / 2, n_keys),
      centers,
      pch = 16,
      col = metrics$colors,
      cex = .cell_relative_size_to_cex(0.5 * key * 25.4),
      xpd = TRUE
    )
    .cell_base_legend_text(
      x = left + key + spacing,
      centers = centers,
      label_lines = metrics$label_lines,
      line_height = label_line,
      cex = metrics$label_cex,
      color = text_color
    )
    return(invisible(NULL))
  }

  bar_width <- geometry$bar_width * unit
  bar_bottom <- body_top - body_height
  n_fills <- length(metrics$fills)
  edges <- seq(bar_bottom, body_top, length.out = n_fills + 1L)
  graphics::rect(
    xleft = left,
    ybottom = edges[-(n_fills + 1L)],
    xright = left + bar_width,
    ytop = edges[-1L],
    col = metrics$fills,
    border = NA,
    xpd = TRUE
  )

  break_y <- scales::rescale(
    metrics$breaks,
    from = metrics$limits,
    to = c(bar_bottom, body_top)
  )
  keep <- .cell_base_legend_label_spacing(
    break_y,
    min_spacing = max_label_lines * label_line
  )
  break_y <- break_y[keep]
  label_lines <- metrics$label_lines[keep]
  tick_length <- geometry$tick_length * bar_width
  graphics::segments(
    x0 = c(rep(left, length(break_y)), rep(left + bar_width - tick_length, length(break_y))),
    y0 = c(break_y, break_y),
    x1 = c(rep(left + tick_length, length(break_y)), rep(left + bar_width, length(break_y))),
    y1 = c(break_y, break_y),
    col = "white",
    lwd = 0.5 * metrics$label_cex,
    xpd = TRUE
  )
  .cell_base_legend_text(
    x = left + bar_width + spacing,
    centers = break_y,
    label_lines = label_lines,
    line_height = label_line,
    cex = metrics$label_cex,
    color = text_color
  )
  return(invisible(NULL))
}

#' Draw wrapped legend labels centered on keys or breaks
#'
#' Each label may occupy several lines. The lines of one label are stacked
#' around its center so a wrapped label stays aligned with its key or tick.
#'
#' @param x Left edge of the label text, in inches.
#' @param centers Vertical center of each label, in inches.
#' @param label_lines A list of character vectors, one vector per label.
#' @param line_height Distance between lines of one label, in inches.
#' @param cex Text size as a `cex` multiplier.
#' @param color Text color.
#'
#' @return `NULL`, invisibly.
#'
#' @noRd
.cell_base_legend_text <- function(
  x,
  centers,
  label_lines,
  line_height,
  cex,
  color
) {
  for (i in seq_along(centers)) {
    lines <- label_lines[[i]]
    n_lines <- length(lines)
    if (n_lines == 0L) {
      next
    }
    ys <- centers[[i]] +
      ((n_lines + 1) / 2 - seq_len(n_lines)) * line_height
    graphics::text(
      x,
      ys,
      labels = lines,
      adj = c(0, 0.5),
      cex = cex,
      col = color,
      xpd = TRUE
    )
  }
  return(invisible(NULL))
}

#' Select legend breaks whose labels do not overlap
#'
#' Walks the break positions from the bottom of the color bar upwards and
#' drops any break closer than `min_spacing` to the previous kept break. The
#' lowest break is always kept.
#'
#' @param positions Increasing numeric break positions in inches.
#' @param min_spacing Minimum distance between labelled breaks in inches.
#'
#' @return A logical vector marking the breaks to label.
#'
#' @noRd
.cell_base_legend_label_spacing <- function(positions, min_spacing) {
  keep <- rep(TRUE, length(positions))
  if (length(positions) < 2L) {
    return(keep)
  }
  last_kept <- positions[[1]]
  for (i in seq_along(positions)[-1L]) {
    if (positions[[i]] - last_kept < min_spacing) {
      keep[[i]] <- FALSE
    } else {
      last_kept <- positions[[i]]
    }
  }
  return(keep)
}

#' Build a continuous base-graphics legend scale
#'
#' Orders colors from low at the bottom to high at the top and places labelled
#' breaks the same way as the ggplot2 color bar: breaks come from
#' [scales::extended_breaks()], breaks outside the scale limits are dropped,
#' and labels use the default ggplot2 number formatting.
#'
#' @param color_scale A built continuous color scale.
#' @param n Number of colors in the legend.
#'
#' @return A list containing `fills`, numeric `breaks`, and character
#' `labels`.
#'
#' @noRd
.cell_base_continuous_legend_scale <- function(color_scale, n = 80L) {
  palette <- .cell_colors_to_hex(unname(color_scale$colors))
  fills <- scales::gradient_n_pal(palette)(seq(0, 1, length.out = n))
  breaks <- .cell_base_legend_breaks(color_scale$limits)
  return(list(
    fills = fills,
    breaks = breaks,
    labels = .cell_base_legend_labels(breaks)
  ))
}

#' Choose legend breaks within continuous color limits
#'
#' @param limits Numeric scale limits of length two.
#'
#' @return Numeric breaks inside `limits`. When no regular break falls inside
#' the limits, the limits themselves are returned.
#'
#' @noRd
.cell_base_legend_breaks <- function(limits) {
  limits <- range(limits)
  breaks <- scales::extended_breaks()(limits)
  breaks <- unique(
    breaks[is.finite(breaks) & breaks >= limits[1] & breaks <= limits[2]]
  )
  if (length(breaks) == 0L) {
    breaks <- unique(limits)
  }
  return(breaks)
}

#' Format legend break labels like ggplot2
#'
#' @param breaks Numeric breaks.
#'
#' @return Character labels with a shared number of decimals.
#'
#' @noRd
.cell_base_legend_labels <- function(breaks) {
  return(format(breaks, trim = TRUE, justify = "left"))
}
