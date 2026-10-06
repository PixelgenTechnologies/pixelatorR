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

  layout <- .cell_base_layout_matrix(
    n_row = n_row,
    n_col = n_col,
    has_title = has_title,
    has_subtitle = has_subtitle,
    has_row_strips = has_row_strips,
    has_col_strips = has_col_strips,
    has_legend = has_legend
  )
  graphics::layout(
    layout$mat,
    widths = layout$widths,
    heights = layout$heights,
    respect = layout$respect
  )

  text_cex <- object$theme$text_size / 12
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
        cex = text_cex,
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
        cex = text_cex,
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
    .cell_base_legend(object, text_cex = text_cex)
  }

  return(invisible(NULL))
}

#' Build a base-graphics layout for a cell plot
#'
#' @param n_row,n_col Panel counts.
#' @param has_title,has_subtitle,has_row_strips,has_col_strips,has_legend
#' Layout flags.
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
  has_legend
) {
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
    if (has_legend) 0.32
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

#' Draw a categorical or continuous legend
#'
#' @param object A `cell_plot_built` object.
#' @param text_cex Text size as a `cex` multiplier.
#'
#' @return `NULL`, invisibly.
#'
#' @noRd
.cell_base_legend <- function(object, text_cex) {
  graphics::par(mar = c(1, 0.4, 1, 0.8))
  graphics::plot.new()
  graphics::rect(
    0,
    0,
    1,
    1,
    col = object$theme$background_color,
    border = NA
  )
  if (identical(object$color$type, "categorical")) {
    graphics::legend(
      "center",
      legend = object$color$limits,
      pch = 16,
      col = unname(object$color$colors[object$color$limits]),
      bty = "n",
      title = .cell_plot_legend_title(object),
      text.col = object$theme$text_color,
      title.col = object$theme$text_color,
      cex = text_cex,
      xpd = TRUE
    )
    return(invisible(NULL))
  }

  legend_scale <- .cell_base_continuous_legend_scale(object$color)
  n <- 80L
  ys <- seq(0.15, 0.85, length.out = n)
  graphics::rect(
    xleft = 0.2,
    ybottom = ys[-n],
    xright = 0.45,
    ytop = ys[-1],
    col = legend_scale$fills[-n],
    border = NA,
    xpd = TRUE
  )
  graphics::text(
    0.55,
    c(0.15, 0.85),
    labels = prettyNum(legend_scale$labels),
    adj = 0,
    cex = text_cex,
    col = object$theme$text_color,
    xpd = TRUE
  )
  graphics::text(
    0.5,
    0.95,
    labels = .cell_plot_legend_title(object),
    cex = text_cex,
    col = object$theme$text_color,
    xpd = TRUE
  )
  return(invisible(NULL))
}

#' Build a continuous base-graphics legend scale
#'
#' Orders colors and endpoint labels from low at the bottom to high at the top,
#' matching the ggplot2 color bar.
#'
#' @param color_scale A built continuous color scale.
#' @param n Number of colors in the legend.
#'
#' @return A list containing `fills` and endpoint `labels`.
#'
#' @noRd
.cell_base_continuous_legend_scale <- function(color_scale, n = 80L) {
  palette <- .cell_colors_to_hex(unname(color_scale$colors))
  fills <- scales::gradient_n_pal(palette)(seq(0, 1, length.out = n))
  return(list(fills = fills, labels = color_scale$limits))
}
