#' Render a cell plot with rgl
#'
#' Builds a [cell_plot()] recipe and draws it as an interactive native 3D
#' scatter using [rgl]. Node sizes are converted from backend-neutral relative
#' units to rgl point diameters in pixels. Continuous sizes are grouped into a
#' bounded number of pixel-size bins to keep the scene responsive.
#' Panel grids use [rgl::layout3d()] with shared mouse control among data
#' panels. Facet labels are drawn in dedicated themeable strip regions along
#' the top (columns) and side (rows), so points cannot cover them. Plot titles
#' and color legends also use reserved regions. Pointer rotation and
#' scroll-zooming stay on the data panels, which move together. Titles, strips,
#' and legends are orthographic and do not rotate or zoom the panels.
#' Panel grids support at most 10 rows and 20 columns.
#'
#' Occlusion follows the scene camera, so markers closer to the current
#' viewpoint appear in front. The rgl backend does not provide hover labels or
#' interactive legend filtering. The `arrange` and `depth` mappings and
#' [cell_coord_rotate()] are ignored. Markers use the rendered colors
#' resolved by [build_cell_plot()], including illumination when requested,
#' while legends are drawn from the unilluminated color scale metadata.
#' Numeric color mappings use a continuous colorbar;
#' categorical mappings use a discrete legend.
#'
#' `output = "window"` opens an rgl device and returns the device id invisibly.
#' `output = "html"` draws the same scene on a null device and returns an rgl
#' htmlwidget for Quarto, R Markdown, and `htmlwidgets::saveWidget()`. Printed
#' at the console, the widget opens in the IDE viewer, fills it, and follows
#' it when the pane is resized. Legends are native scene objects and do not
#' support Plotly-style interactive legend filtering. Title and legend text
#' uses rgl's own font. Facet strip labels are rasterized with base graphics,
#' because rgl text cannot be rotated for the row strips.
#'
#' @param object A `cell_plot` recipe.
#' @param output `"window"` opens a local rgl window. `"html"` returns an
#' htmlwidget.
#' @param width,height Size in pixels. For `"window"` this is the window
#' size, 1000 by 1000 when not given. For `"html"` it is the canvas size.
#' Inside a knitr HTML chunk, a missing html size is the chunk `fig.width` or
#' `fig.height` in inches multiplied by `dpi`. Outside knitr, an html widget
#' without `width` and `height` fills the viewer or browser element it is
#' shown in, and its layout is computed for a 1000 by 1000 canvas.
#'
#' @return For `"window"`, the rgl device id, invisibly. For `"html"`, an
#' htmlwidget.
#'
#' @seealso [cell_plot()], [cell_plot_interactive()]
#'
#' @examplesIf interactive()
#' se <- ReadPNA_Seurat(minimal_pna_pxl_file())
#' se <- LoadCellGraphs(se, cells = colnames(se)[4], verbose = FALSE) |>
#'   ComputeLayout(layout_method = "spectral")
#'
#' cell_graph <- CellGraphs(se)[[4]]
#'
#' layout_data <- FetchLayoutData(cell_graph, vars = "CD82", layout_method = "spectral_3d")
#'
#' cell_plot(layout_data, color = CD82) |>
#'   cell_plot_rgl()
#'
#' @export
cell_plot_rgl <- function(
  object,
  output = c("window", "html"),
  width = NULL,
  height = NULL
) {
  .validate_cell_plot(object)
  expect_rgl()
  output <- rlang::arg_match(output)
  width <- .cell_rgl_check_px(width, arg = "width")
  height <- .cell_rgl_check_px(height, arg = "height")

  object$mapping$arrange <- NULL
  built <- build_cell_plot(object)
  if (identical(output, "html")) {
    return(.cell_rgl_html_widget(built, width = width, height = height))
  }
  size <- .cell_rgl_canvas_pixels(width, height)
  return(.render_cell_plot_rgl(built, size = size))
}

#' Draw a cell plot into an rgl htmlwidget
#'
#' Opens a null device, draws the built plot, snapshots it with
#' [rgl::rglwidget()], and closes the device. The snapshot keeps the scene
#' after the device closes. When neither the caller nor a knitr HTML chunk
#' sets a size, the widget carries no size of its own, so htmlwidgets lets it
#' fill the IDE viewer and follow the viewer when it is resized.
#'
#' @param object A `cell_plot_built` object.
#' @param width,height Checked canvas sizes in pixels, or `NULL`.
#'
#' @return An rgl htmlwidget.
#'
#' @noRd
.cell_rgl_html_widget <- function(object, width, height) {
  chunk <- .cell_rgl_knitr_html_px()
  sized <- !is.null(width) || !is.null(height) || !is.null(chunk)
  size <- .cell_rgl_canvas_pixels(width, height, fallback = chunk)
  previous <- as.integer(rgl::cur3d())
  device <- rgl::open3d(
    useNULL = TRUE,
    silent = TRUE,
    windowRect = c(0, 0, size$width, size$height)
  )
  on.exit(.cell_rgl_close_html_device(device, previous), add = TRUE)
  .render_cell_plot_rgl(object, device = device, webgl = TRUE)
  # A knitr PDF or Word chunk would otherwise ask for a raster snapshot.
  # A null device cannot draw one, and this output is the widget itself.
  return(rgl::rglwidget(
    width = if (sized) size$width,
    height = if (sized) size$height,
    snapshot = FALSE
  ))
}

#' Close the null device used for an rgl htmlwidget
#'
#' Restores the device that was current before the widget was drawn, when that
#' device is still open.
#'
#' @param device Device id opened for the widget.
#' @param previous Device id that was current beforehand.
#'
#' @return `NULL`, invisibly.
#'
#' @noRd
.cell_rgl_close_html_device <- function(device, previous) {
  open_devices <- as.integer(rgl::rgl.dev.list())
  if (as.integer(device) %in% open_devices) {
    rgl::set3d(device, silent = TRUE)
    rgl::close3d()
  }
  open_devices <- as.integer(rgl::rgl.dev.list())
  if (length(previous) == 1L && previous %in% open_devices) {
    rgl::set3d(previous, silent = TRUE)
  }
  return(invisible(NULL))
}

#' Pixel size of an rgl canvas
#'
#' An explicit size wins, then the `fallback`, then a 1000 by 1000 square.
#'
#' @param width,height Checked sizes in pixels, or `NULL`.
#' @param fallback Optional list with `width` and `height` used for a missing
#' dimension, such as the knitr chunk figure size.
#'
#' @return A list with integer `width` and `height`.
#'
#' @noRd
.cell_rgl_canvas_pixels <- function(width, height, fallback = NULL) {
  if (is.null(width)) {
    width <- fallback$width %||% .cell_rgl_default_canvas_px
  }
  if (is.null(height)) {
    height <- fallback$height %||% .cell_rgl_default_canvas_px
  }
  return(list(width = as.integer(width), height = as.integer(height)))
}

#' Default rgl canvas size in pixels
#'
#' Used for both the local window and the html canvas when no size is given.
#'
#' @noRd
.cell_rgl_default_canvas_px <- 1000L

#' Check one canvas dimension
#'
#' @param value A pixel count, or `NULL` when the caller did not set it.
#' @param arg Argument name used in the error.
#' @param call Calling environment used for validation errors.
#'
#' @return `NULL`, or one positive integer pixel count.
#'
#' @noRd
.cell_rgl_check_px <- function(value, arg, call = rlang::caller_env()) {
  if (is.null(value)) {
    return(NULL)
  }
  assert_single_value(value, type = "integer", arg = arg, call = call)
  if (!is.finite(value) || value <= 0) {
    cli::cli_abort(
      c("x" = "{.arg {arg}} must be a positive number of pixels."),
      call = call
    )
  }
  return(as.integer(value))
}

#' Figure size of the current knitr HTML chunk
#'
#' Returns pixel sizes only while knitr is rendering an HTML document. Other
#' outputs, including PDF and Word, keep the default canvas.
#'
#' @return A list with integer `width` and `height`, or `NULL` when the chunk
#' size does not apply.
#'
#' @noRd
.cell_rgl_knitr_html_px <- function() {
  if (!isTRUE(getOption("knitr.in.progress"))) {
    return(NULL)
  }
  if (!requireNamespace("knitr", quietly = TRUE)) {
    return(NULL)
  }
  if (!isTRUE(knitr::is_html_output())) {
    return(NULL)
  }
  fig_width <- knitr::opts_current$get("fig.width")
  fig_height <- knitr::opts_current$get("fig.height")
  dpi <- knitr::opts_current$get("dpi")
  if (is.null(fig_width) || is.null(fig_height) || is.null(dpi)) {
    return(NULL)
  }
  return(list(
    width = as.integer(round(fig_width * dpi)),
    height = as.integer(round(fig_height * dpi))
  ))
}

#' Render a built cell plot with rgl
#'
#' Converts rendered colors, sizes, and alpha into rgl point properties.
#' Colors come from .cell_plot_rendered_colors() so builder-baked illumination
#' matches ggplot, Plotly, and base. Panel grids become an [rgl::layout3d()]
#' arrangement that keeps every row and column level, including empty panels.
#' Facet labels, titles, and legends are native objects in orthographic
#' subscenes, so placement does not depend on the upper-left panel.
#'
#' @param object A `cell_plot_built` object with a mapped `z` coordinate.
#' @param device An open rgl device to draw into. When `NULL`, a local window
#' is opened.
#' @param size List with `width` and `height` in pixels for a new window.
#' Ignored when `device` is given.
#' @param webgl Whether the scene will be shown through [rgl::rglwidget()]
#' rather than in the window. Chrome camera fitting differs between the two.
#'
#' @return The rgl device id, returned invisibly.
#'
#' @noRd
.render_cell_plot_rgl <- function(
  object,
  device = NULL,
  size = .cell_rgl_canvas_pixels(NULL, NULL),
  webgl = FALSE
) {
  pixelatorR:::assert_class(object, "cell_plot_built", arg = "object")

  mapping <- object$mapping
  categorical <- identical(object$color$type, "categorical")
  n_nodes <- nrow(object$data)

  plot_data <- object$data
  # Illumination is pre-baked into color$illuminated by build_cell_plot();
  # use the same rendered colors as ggplot/Plotly/base, not the raw scale.
  plot_data$.cell_color <- rep_len(.cell_plot_rendered_colors(object), n_nodes)
  plot_data$.cell_alpha <- rep_len(object$alpha$resolved, n_nodes)
  # Depth sizing is ggplot-only; rgl uses resolved sizes as-is.
  plot_data$.cell_size <- .cell_relative_size_to_pixels(
    rep_len(object$size$resolved, n_nodes)
  )

  ranges <- list(
    x = .cell_plot_padded_range(plot_data[[mapping$x]]),
    y = .cell_plot_padded_range(plot_data[[mapping$y]]),
    z = .cell_plot_padded_range(plot_data[[mapping$z]])
  )

  facet_rows <- object$grid$rows
  facet_cols <- object$grid$cols
  row_levels <- if (is.null(facet_rows)) {
    list("")
  } else {
    .cell_plot_facet_levels(plot_data[[facet_rows]])
  }
  col_levels <- if (is.null(facet_cols)) {
    list("")
  } else {
    .cell_plot_facet_levels(plot_data[[facet_cols]])
  }
  if (length(row_levels) > 10) {
    cli::cli_abort(
      c(
        "x" = "{.fn cell_plot_rgl} supports at most 10 facet rows.",
        "i" = "{.arg rows} creates {length(row_levels)} facet rows."
      )
    )
  }
  if (length(col_levels) > 20) {
    cli::cli_abort(
      c(
        "x" = "{.fn cell_plot_rgl} supports at most 20 facet columns.",
        "i" = "{.arg cols} creates {length(col_levels)} facet columns."
      )
    )
  }

  n_row <- length(row_levels)
  n_col <- length(col_levels)
  n_panels <- n_row * n_col
  background_color <- object$theme$background_color
  text_color <- object$theme$text_color
  text_size <- object$theme$text_size
  strip_background_color <- object$theme$strip_background_color

  title <- object$annotation$title
  subtitle <- object$annotation$subtitle
  need_title <- !is.null(title) || !is.null(subtitle)
  need_legend <- !is.null(mapping$color)
  has_layout <- !is.null(facet_rows) ||
    !is.null(facet_cols) ||
    need_title ||
    need_legend
  legend_title <- .cell_plot_legend_title(object)
  categorical_legend <- if (need_legend && categorical) {
    list(
      title = legend_title,
      labels = object$color$limits,
      colors = .cell_colors_to_hex(
        unname(object$color$colors[object$color$limits])
      )
    )
  }
  continuous_legend <- if (need_legend && !categorical) {
    list(
      title = legend_title,
      limits = object$color$limits,
      colors = .cell_colors_to_hex(unname(object$color$colors))
    )
  }

  if (is.null(device)) {
    device <- rgl::open3d(
      windowRect = c(100, 100, 100 + size$width, 100 + size$height)
    )
  } else {
    rgl::set3d(device, silent = TRUE)
  }

  title_id <- NULL
  legend_id <- NULL
  corner_id <- NULL
  col_strip_ids <- integer()
  row_strip_ids <- integer()
  if (!has_layout) {
    panel_ids <- rgl::subsceneInfo()$id
  } else {
    parent_id <- rgl::currentSubscene3d()
    layout <- do.call(
      .cell_rgl_facet_layout,
      list(
        n_row = n_row,
        n_col = n_col,
        need_col_strips = !is.null(facet_cols),
        need_row_strips = !is.null(facet_rows),
        need_title = need_title,
        need_legend = need_legend,
        need_subtitle = !is.null(subtitle),
        viewport = rgl::par3d("viewport", subscene = parent_id)
      )
    )
    # mouseMode = "replace" is required: layout3d() defaults to inherited
    # mouse handling, so disabling chrome mouse would write through to the
    # parent and kill trackball/zoom on every data panel.
    ids <- rgl::layout3d(
      layout$mat,
      widths = layout$widths,
      heights = layout$heights,
      sharedMouse = FALSE,
      mouseMode = "replace"
    )
    panel_ids <- ids[seq_len(n_panels)]
    chrome_index <- n_panels
    if (!is.null(facet_cols)) {
      col_strip_ids <- ids[chrome_index + seq_len(n_col)]
      chrome_index <- chrome_index + n_col
    }
    if (!is.null(facet_cols) && !is.null(facet_rows)) {
      chrome_index <- chrome_index + 1L
      corner_id <- ids[[chrome_index]]
    }
    if (!is.null(facet_rows)) {
      row_strip_ids <- ids[chrome_index + seq_len(n_row)]
      chrome_index <- chrome_index + n_row
    }
    if (need_title) {
      chrome_index <- chrome_index + 1L
      title_id <- ids[[chrome_index]]
    }
    if (need_legend) {
      chrome_index <- chrome_index + 1L
      legend_id <- ids[[chrome_index]]
    }
    .cell_rgl_set_listeners(panel_ids, panel_ids = panel_ids)
  }

  for (panel_index in seq_len(n_panels)) {
    # Select panels by id so title/legend chrome cells in the layout3d
    # subscene list cannot shift which data panel is current.
    rgl::useSubscene3d(panel_ids[[panel_index]])

    row_index <- ((panel_index - 1L) %/% n_col) + 1L
    col_index <- ((panel_index - 1L) %% n_col) + 1L

    keep <- rep(TRUE, n_nodes)
    if (!is.null(facet_rows)) {
      keep <- keep &
        .cell_plot_facet_match(
          plot_data[[facet_rows]],
          row_levels[[row_index]]
        )
    }
    if (!is.null(facet_cols)) {
      keep <- keep &
        .cell_plot_facet_match(
          plot_data[[facet_cols]],
          col_levels[[col_index]]
        )
    }
    panel_data <- plot_data[keep, , drop = FALSE]

    # Establish a shared data-aspect bounding box for every panel, including
    # empty ones, matching facet_grid(drop = FALSE) and Plotly scenes.
    rgl::plot3d(
      x = ranges$x,
      y = ranges$y,
      z = ranges$z,
      type = "n",
      xlab = "",
      ylab = "",
      zlab = "",
      axes = FALSE,
      decorate = FALSE
    )

    if (nrow(panel_data) > 0) {
      .cell_rgl_points(
        x = panel_data[[mapping$x]],
        y = panel_data[[mapping$y]],
        z = panel_data[[mapping$z]],
        color = panel_data$.cell_color,
        size = panel_data$.cell_size,
        alpha = panel_data$.cell_alpha
      )
    }

    rgl::bg3d(color = background_color)
    rgl::view3d(zoom = 1)
  }

  for (col_index in seq_along(col_strip_ids)) {
    .cell_rgl_draw_strip(
      subscene = col_strip_ids[[col_index]],
      label = .cell_plot_facet_label(col_levels[[col_index]]),
      angle = 0,
      text_color = text_color,
      text_size = text_size,
      background_color = strip_background_color,
      webgl = webgl
    )
  }
  for (row_index in seq_along(row_strip_ids)) {
    .cell_rgl_draw_strip(
      subscene = row_strip_ids[[row_index]],
      label = .cell_plot_facet_label(row_levels[[row_index]]),
      angle = .cell_plot_row_strip_angle,
      text_color = text_color,
      text_size = text_size,
      background_color = strip_background_color,
      webgl = webgl
    )
  }
  if (!is.null(corner_id)) {
    .cell_rgl_draw_strip(
      subscene = corner_id,
      label = "",
      angle = 0,
      text_color = text_color,
      text_size = text_size,
      background_color = background_color,
      webgl = webgl
    )
  }
  if (!is.null(title_id)) {
    .cell_rgl_draw_title(
      subscene = title_id,
      title = title,
      subtitle = subtitle,
      text_color = text_color,
      text_size = text_size,
      background_color = background_color,
      webgl = webgl
    )
  }
  if (!is.null(legend_id)) {
    .cell_rgl_draw_legend(
      subscene = legend_id,
      text_color = text_color,
      text_size = text_size,
      background_color = background_color,
      categorical_legend = categorical_legend,
      continuous_legend = continuous_legend,
      webgl = webgl
    )
  }

  return(invisible(as.integer(device)))
}

#' Build a layout matrix for faceted rgl chrome
#'
#' Reserves dedicated column and row facet strips around an `n_row` by `n_col`
#' panel grid, plus optional title and legend regions. When both strips are
#' present, the strip-corner cell is also reserved so pointer events there
#' reach a subscene instead of the unwired root. Panel cells are numbered
#' first, followed by column strips, the strip corner, row strips, title, and
#' legend, so [rgl::layout3d()] returns ids in that order.
#'
#' @param n_row,n_col Panel grid dimensions.
#' @param need_col_strips,need_row_strips Whether to reserve facet strips.
#' @param need_title,need_legend Whether to reserve title and legend regions.
#' @param need_subtitle Whether the title region also contains a subtitle.
#' @param viewport Parent viewport as either width and height or the four-value
#' rgl viewport vector.
#'
#' @return A list with `mat`, `widths`, and `heights`.
#'
#' @noRd
.cell_rgl_facet_layout <- function(
  n_row,
  n_col,
  need_col_strips,
  need_row_strips,
  need_title,
  need_legend,
  need_subtitle = FALSE,
  viewport = c(width = 1, height = 1)
) {
  n_panels <- n_row * n_col
  layout_rows <- n_row +
    as.integer(need_col_strips) +
    as.integer(need_title)
  layout_cols <- n_col +
    as.integer(need_row_strips) +
    as.integer(need_legend)
  mat <- matrix(0L, nrow = layout_rows, ncol = layout_cols)

  panel_index <- 0L
  row_offset <- as.integer(need_title) + as.integer(need_col_strips)
  col_offset <- as.integer(need_row_strips)
  for (row_index in seq_len(n_row)) {
    for (col_index in seq_len(n_col)) {
      panel_index <- panel_index + 1L
      mat[row_index + row_offset, col_index + col_offset] <- panel_index
    }
  }

  chrome_id <- n_panels
  if (need_col_strips) {
    strip_row <- as.integer(need_title) + 1L
    for (col_index in seq_len(n_col)) {
      chrome_id <- chrome_id + 1L
      mat[strip_row, col_index + col_offset] <- chrome_id
    }
  }
  if (need_col_strips && need_row_strips) {
    chrome_id <- chrome_id + 1L
    mat[as.integer(need_title) + 1L, 1L] <- chrome_id
  }
  if (need_row_strips) {
    for (row_index in seq_len(n_row)) {
      chrome_id <- chrome_id + 1L
      mat[row_index + row_offset, 1L] <- chrome_id
    }
  }
  if (need_title) {
    chrome_id <- chrome_id + 1L
    mat[1L, seq_len(n_col + col_offset)] <- chrome_id
  }
  if (need_legend) {
    chrome_id <- chrome_id + 1L
    mat[, layout_cols] <- chrome_id
  }

  widths <- c(
    if (need_row_strips) 0.1,
    rep(1, n_col),
    if (need_legend) 0.28
  )
  heights <- c(
    if (need_title) {
      if (need_subtitle) 0.18 else 0.12
    },
    if (need_col_strips) 0.1,
    rep(1, n_row)
  )
  if (need_col_strips && need_row_strips) {
    viewport <- as.numeric(viewport)
    if (length(viewport) == 4L) {
      viewport <- viewport[3:4]
    }
    viewport_width <- viewport[[1]]
    viewport_height <- viewport[[2]]
    col_strip_row <- as.integer(need_title) + 1L
    col_strip_height <- viewport_height *
      heights[[col_strip_row]] / sum(heights)
    other_width_weight <- sum(widths[-1L])
    if (col_strip_height < viewport_width) {
      widths[[1L]] <- col_strip_height * other_width_weight /
        (viewport_width - col_strip_height)
    }
  }
  return(list(mat = mat, widths = widths, heights = heights))
}

#' Connect rgl data panels so they share one camera
#'
#' Each data panel listens to every data panel. Title, strip, and legend
#' subscenes are left alone, so pointer movement there does not move the plot.
#'
#' @param subscene_ids Integer vector of subscene ids to configure.
#' @param panel_ids Integer vector of data-panel subscene ids.
#'
#' @return `subscene_ids`, invisibly.
#'
#' @noRd
.cell_rgl_set_listeners <- function(subscene_ids, panel_ids) {
  for (id in subscene_ids) {
    rgl::par3d(listeners = panel_ids, subscene = id)
  }
  return(invisible(subscene_ids))
}

#' Prepare an orthographic rgl chrome subscene
#'
#' Titles, strips, and legends keep a fixed camera. `FOV = 0` is orthographic,
#' and every mouse button is `none`, so the region cannot rotate or zoom the
#' data panels.
#'
#' Chrome is laid out in the region's own pixel coordinates: `x` runs from 0
#' to the viewport width and `y` from 0 to its height. A zero-length segment
#' pair at two opposite corners sets that bounding box without drawing a
#' visible frame. rgl fits the bounding sphere of the box to the shorter
#' viewport side, which leaves a tall or wide region mostly empty, so the
#' zoom is set to the factor that makes the box itself fill the viewport.
#' WebGL encloses a sphere 1.1 times larger than the window does, which the
#' `webgl` flag compensates for.
#'
#' @param subscene Subscene id to draw in.
#' @param background_color Background color for the region.
#' @param webgl Whether the scene is drawn for an htmlwidget.
#'
#' @return A list with the viewport `width` and `height` in pixels.
#'
#' @noRd
.cell_rgl_prepare_chrome <- function(subscene, background_color, webgl) {
  rgl::useSubscene3d(subscene)
  viewport <- as.numeric(rgl::par3d("viewport", subscene = subscene))
  width <- max(viewport[[3]], 1)
  height <- max(viewport[[4]], 1)
  rgl::bg3d(color = background_color)
  rgl::par3d(
    FOV = 0,
    mouseMode = c(
      left = "none",
      right = "none",
      middle = "none",
      wheel = "none"
    ),
    subscene = subscene
  )
  rgl::segments3d(
    x = c(0, 0, width, width),
    y = c(0, 0, height, height),
    z = c(0, 0, 0, 0),
    color = background_color,
    lit = FALSE
  )
  # The window maps the half-diagonal of the box to the shorter half-side of
  # the viewport. Zooming by shorter side over diagonal maps the box edge
  # there instead.
  zoom <- min(width, height) / sqrt(width^2 + height^2)
  if (webgl) {
    zoom <- zoom / .cell_rgl_webgl_enclose_factor
  }
  rgl::view3d(theta = 0, phi = 0, fov = 0, zoom = zoom)
  return(list(width = width, height = height))
}

#' Orthographic fit ratio between rgl's WebGL and window renderers
#'
#' rgl's WebGL projection encloses a bounding sphere 1.1 times the scene
#' radius while keeping the same observer distance, so its orthographic
#' half-length is 0.9 times the one the window uses.
#'
#' @noRd
.cell_rgl_webgl_enclose_factor <- 0.9

#' Convert a theme text size into an rgl cex
#'
#' Theme text size uses the same point size as the ggplot theme, where 11 is
#' the default. rgl draws its own font at `cex = 1` for that default.
#'
#' @param text_size Theme text size.
#'
#' @return A cex multiplier.
#'
#' @noRd
.cell_rgl_text_cex <- function(text_size) {
  return(text_size / 11)
}

#' Draw a facet label in a dedicated rgl strip
#'
#' The strip is a separate subscene from the data panel, so points cannot cover
#' its text. A gray background follows the default [ggplot2::facet_grid()]
#' appearance. rgl text always faces the screen upright, so labels are drawn
#' through [rgl::plotmath3d()], which rasterizes the label with base graphics
#' at the requested angle and places it as a fixed-size sprite. Column strips
#' use the same path so both strip kinds share one font.
#'
#' @param subscene Subscene id for the strip.
#' @param label Facet level label.
#' @param angle Text rotation in degrees.
#' @param text_color,text_size,background_color Theme values.
#' @param webgl Whether the scene is drawn for an htmlwidget.
#'
#' @return `NULL`, invisibly.
#'
#' @noRd
.cell_rgl_draw_strip <- function(
  subscene,
  label,
  angle,
  text_color,
  text_size,
  background_color,
  webgl
) {
  region <- .cell_rgl_prepare_chrome(
    subscene = subscene,
    background_color = background_color,
    webgl = webgl
  )
  if (!nzchar(label)) {
    return(invisible(NULL))
  }
  rgl::plotmath3d(
    x = region$width / 2,
    y = region$height / 2,
    z = 0,
    text = label,
    adj = 0.5,
    cex = .cell_rgl_text_cex(text_size),
    srt = angle,
    col = text_color
  )
  return(invisible(NULL))
}

#' Draw a title in a reserved rgl layout region
#'
#' The title and subtitle are left-aligned at the region's left edge, as in
#' ggplot, with a margin of one text size in pixels.
#'
#' @param subscene Subscene id for the title.
#' @param title,subtitle Optional plot title and subtitle text.
#' @param text_color,text_size,background_color Theme values.
#' @param webgl Whether the scene is drawn for an htmlwidget.
#'
#' @return `NULL`, invisibly.
#'
#' @noRd
.cell_rgl_draw_title <- function(
  subscene,
  title,
  subtitle,
  text_color,
  text_size,
  background_color,
  webgl
) {
  region <- .cell_rgl_prepare_chrome(
    subscene = subscene,
    background_color = background_color,
    webgl = webgl
  )
  cex <- .cell_rgl_text_cex(text_size)
  left <- text_size
  if (!is.null(title)) {
    rgl::text3d(
      x = left,
      y = region$height * if (is.null(subtitle)) 0.5 else 0.64,
      z = 0,
      texts = title,
      adj = c(0, 0.5),
      color = text_color,
      cex = cex * 1.2,
      font = 2
    )
  }
  if (!is.null(subtitle)) {
    rgl::text3d(
      x = left,
      y = region$height * if (is.null(title)) 0.5 else 0.28,
      z = 0,
      texts = subtitle,
      adj = c(0, 0.5),
      color = text_color,
      cex = cex
    )
  }
  return(invisible(NULL))
}

#' Draw a legend in a reserved rgl layout region
#'
#' @param subscene Subscene id for the legend.
#' @param text_color,text_size,background_color Theme values.
#' @param categorical_legend Optional list with `title`, `labels`, and `colors`.
#' @param continuous_legend Optional list with `title`, `limits`, and `colors`.
#' @param webgl Whether the scene is drawn for an htmlwidget.
#'
#' @return `NULL`, invisibly.
#'
#' @noRd
.cell_rgl_draw_legend <- function(
  subscene,
  text_color,
  text_size,
  background_color,
  categorical_legend = NULL,
  continuous_legend = NULL,
  webgl
) {
  region <- .cell_rgl_prepare_chrome(
    subscene = subscene,
    background_color = background_color,
    webgl = webgl
  )
  if (!is.null(continuous_legend)) {
    .cell_rgl_draw_colorbar(
      continuous_legend = continuous_legend,
      text_color = text_color,
      text_size = text_size,
      region = region
    )
  } else if (!is.null(categorical_legend)) {
    .cell_rgl_draw_discrete_legend(
      categorical_legend = categorical_legend,
      text_color = text_color,
      text_size = text_size,
      region = region
    )
  }
  return(invisible(NULL))
}

#' Draw a categorical legend with point swatches
#'
#' Rows are spaced by the text size in pixels and the block, including its
#' title, is centered vertically in the legend region.
#'
#' @param categorical_legend List with `title`, `labels`, and `colors`.
#' @param text_color Theme text color.
#' @param text_size Theme text size.
#' @param region List with the region `width` and `height` in pixels.
#'
#' @return `NULL`, invisibly.
#'
#' @noRd
.cell_rgl_draw_discrete_legend <- function(
  categorical_legend,
  text_color,
  text_size,
  region
) {
  labels <- categorical_legend$labels
  n_labels <- length(labels)
  if (n_labels == 0L) {
    return(invisible(NULL))
  }
  title <- categorical_legend$title
  has_title <- !is.null(title) && nzchar(title)
  line <- text_size * 1.8
  n_rows <- n_labels + if (has_title) 1.4 else 0
  top <- region$height / 2 + n_rows * line / 2
  swatch_x <- region$width * 0.16
  label_x <- region$width * 0.28
  if (has_title) {
    rgl::text3d(
      x = swatch_x,
      y = top - line * 0.5,
      z = 0,
      texts = title,
      adj = c(0, 0.5),
      color = text_color,
      cex = .cell_rgl_text_cex(text_size)
    )
    top <- top - line * 1.4
  }
  ys <- top - (seq_len(n_labels) - 0.5) * line
  rgl::points3d(
    x = rep(swatch_x, n_labels),
    y = ys,
    z = rep(0, n_labels),
    color = categorical_legend$colors,
    size = text_size * 0.8,
    lit = FALSE
  )
  rgl::text3d(
    x = rep(label_x, n_labels),
    y = ys,
    z = rep(0, n_labels),
    texts = labels,
    adj = c(0, 0.5),
    color = text_color,
    cex = .cell_rgl_text_cex(text_size)
  )
  return(invisible(NULL))
}

#' Limits used to place a continuous colorbar
#'
#' Constant or non-finite limits cannot position ticks. They expand to a short
#' range around the repeated value, or to `c(0, 1)` when nothing is finite.
#'
#' @param limits Numeric color-scale limits.
#'
#' @return A length-two finite, strictly increasing range.
#'
#' @noRd
.cell_rgl_colorbar_limits <- function(limits) {
  if (length(limits) < 2L || any(!is.finite(limits))) {
    return(c(0, 1))
  }
  if (limits[[1]] == limits[[2]]) {
    pad <- max(abs(limits[[1]]) * 0.05, 1e-6)
    return(c(limits[[1]] - pad, limits[[2]] + pad))
  }
  return(limits)
}

#' Draw a continuous colorbar for a numeric color scale
#'
#' Numeric mappings render as stacked [rgl::quads3d()] with vertex colors,
#' rather than a discrete key. Ticks are [rgl::segments3d()] and labels are
#' [rgl::text3d()]. Positions are fractions of the legend region's width and
#' height in pixels.
#'
#' @param continuous_legend List with `title`, `limits`, and `colors`.
#' @param text_color Theme text color.
#' @param text_size Theme text size.
#' @param region List with the region `width` and `height` in pixels.
#'
#' @return `NULL`, invisibly.
#'
#' @noRd
.cell_rgl_draw_colorbar <- function(
  continuous_legend,
  text_color,
  text_size,
  region
) {
  limits <- .cell_rgl_colorbar_limits(continuous_legend$limits)
  colors <- continuous_legend$colors
  if (length(colors) == 0L) {
    colors <- "#000000"
  }
  bar_x <- c(0.18, 0.3) * region$width
  bar_y <- c(0.3, 0.7) * region$height
  if (length(colors) == 1L) {
    rgl::quads3d(
      x = c(bar_x[[1]], bar_x[[2]], bar_x[[2]], bar_x[[1]]),
      y = c(bar_y[[1]], bar_y[[1]], bar_y[[2]], bar_y[[2]]),
      z = rep(0, 4),
      color = colors,
      lit = FALSE
    )
  } else {
    ys <- seq(bar_y[[1]], bar_y[[2]], length.out = length(colors))
    for (stop_index in seq_len(length(colors) - 1L)) {
      rgl::quads3d(
        x = c(bar_x[[1]], bar_x[[2]], bar_x[[2]], bar_x[[1]]),
        y = c(
          ys[[stop_index]],
          ys[[stop_index]],
          ys[[stop_index + 1L]],
          ys[[stop_index + 1L]]
        ),
        z = rep(0, 4),
        color = c(
          colors[[stop_index]],
          colors[[stop_index]],
          colors[[stop_index + 1L]],
          colors[[stop_index + 1L]]
        ),
        lit = FALSE
      )
    }
  }

  tick_values <- pretty(limits, n = 4)
  tick_values <- tick_values[
    tick_values >= limits[[1]] & tick_values <= limits[[2]]
  ]
  if (length(tick_values) > 0L) {
    tick_y <- stats::approx(x = limits, y = bar_y, xout = tick_values)$y
    tick_length <- diff(bar_x) * 0.2
    n_ticks <- length(tick_values)
    rgl::segments3d(
      x = as.vector(rbind(
        rep(bar_x[[2]], n_ticks),
        rep(bar_x[[2]] + tick_length, n_ticks)
      )),
      y = as.vector(rbind(tick_y, tick_y)),
      z = rep(0, 2 * n_ticks),
      color = text_color,
      lit = FALSE
    )
    rgl::text3d(
      x = rep(bar_x[[2]] + tick_length * 1.4, n_ticks),
      y = tick_y,
      z = rep(0, n_ticks),
      texts = format(tick_values, trim = TRUE, scientific = FALSE),
      adj = c(0, 0.5),
      color = text_color,
      cex = .cell_rgl_text_cex(text_size) * 0.85
    )
  }

  title <- continuous_legend$title
  if (!is.null(title) && nzchar(title)) {
    rgl::text3d(
      x = mean(bar_x),
      y = 0.77 * region$height,
      z = 0,
      texts = title,
      color = text_color,
      cex = .cell_rgl_text_cex(text_size) * 0.9
    )
  }
  return(invisible(NULL))
}

#' Maximum number of rgl point-size groups
#'
#' Bounds the number of scene objects produced by a continuous size mapping.
#'
#' @noRd
.cell_rgl_max_size_groups <- 20L

#' Group rgl point sizes into a bounded number of draw calls
#'
#' Exact sizes are retained when there are few distinct values. Continuous
#' mappings with many distinct values are quantized into equal-width bins,
#' preventing one rgl scene object per node while preserving the mapped size
#' range. Missing sizes form one additional fallback group.
#'
#' @param size Numeric rgl point sizes.
#' @param max_groups Maximum number of finite size groups.
#'
#' @return A list containing a group key per point and one representative size
#'   per group.
#'
#' @noRd
.cell_rgl_size_groups <- function(
  size,
  max_groups = .cell_rgl_max_size_groups
) {
  finite_size <- size[is.finite(size)]
  unique_size <- unique(finite_size)
  if (length(unique_size) <= max_groups) {
    keys <- ifelse(is.finite(size), match(size, unique_size), 0L)
  } else {
    size_range <- range(finite_size)
    scaled <- (size[is.finite(size)] - size_range[[1]]) / diff(size_range)
    bin <- as.integer(floor(scaled * max_groups)) + 1L
    bin[bin < 1L] <- 1L
    bin[bin > max_groups] <- max_groups
    keys <- rep(0L, length(size))
    keys[is.finite(size)] <- bin
  }
  representative <- vapply(
    unique(keys),
    function(key) {
      if (!is.finite(key) || key == 0L) {
        return(1)
      }
      return(stats::median(size[keys == key], na.rm = TRUE))
    },
    numeric(1)
  )
  names(representative) <- as.character(unique(keys))
  return(list(keys = keys, size = representative))
}

#' Draw rgl points grouped by pixel size
#'
#' Uses the same unlit `points` primitive as [rgl::plot3d()] `type = "p"`.
#' [rgl::points3d()] accepts only a scalar point size, so nodes are split into
#' a bounded number of size groups while
#' preserving per-node color and alpha vectors within each group.
#'
#' @param x,y,z Node coordinates.
#' @param color Hex colors.
#' @param size rgl point sizes in pixels.
#' @param alpha Per-node opacities.
#'
#' @return Invisibly, the object ids created by [rgl::points3d()].
#'
#' @noRd
.cell_rgl_points <- function(x, y, z, color, size, alpha) {
  ids <- integer()
  groups <- .cell_rgl_size_groups(size)
  for (key in unique(groups$keys)) {
    idx <- which(groups$keys == key)
    ids <- c(
      ids,
      rgl::points3d(
        x = x[idx],
        y = y[idx],
        z = z[idx],
        size = groups$size[[as.character(key)]],
        color = color[idx],
        alpha = alpha[idx]
      )
    )
  }
  return(invisible(ids))
}
