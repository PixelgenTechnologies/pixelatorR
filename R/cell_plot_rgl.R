#' Render a cell plot with rgl
#'
#' `r lifecycle::badge("experimental")`
#'
#' Builds a [cell_plot()] recipe and draws it as an interactive native 3D
#' scatter using [rgl]. Node sizes are converted from backend-neutral relative
#' units to rgl point diameters in pixels. Continuous sizes are grouped into a
#' bounded number of pixel-size bins to keep the scene responsive.
#' Panel grids use [rgl::layout3d()] with shared mouse control among data
#' panels. Facet labels are drawn in dedicated themeable strip regions along
#' the top (columns) and side (rows), so points cannot cover them. Plot titles
#' and color legends also use reserved full-window chrome regions. Rotating and
#' scroll-zooming work with the pointer anywhere in the window, including over
#' the strips, title, and legend, and always drive every data panel together.
#' Panel grids support at most 10 rows and 20 columns.
#'
#' Occlusion follows the scene camera, so markers closer to the current
#' viewpoint appear in front. The rgl backend does not provide hover labels or
#' interactive legend filtering. The `arrange` and `depth` mappings and
#' [cell_coord_rotate()] are ignored. Markers use the rendered colors
#' resolved by [build_cell_plot()], including illumination when requested,
#' while legends are drawn from the unilluminated color scale metadata as
#' static 2D overlays. Numeric color mappings use a continuous colorbar;
#' categorical mappings use a discrete legend.
#'
#' Unlike [cell_plot_interactive()], this renderer opens an rgl device rather
#' than returning an htmlwidget. Legends are drawn with [rgl::bgplot3d()] and
#' do not support Plotly-style interactive legend filtering.
#'
#' [rgl::bgplot3d()] renders text into a bitmap sized for the window it was
#' drawn in. Titles, facet strips, and legends are redrawn automatically at
#' the new size shortly after the window stops changing.
#'
#' @param object A `cell_plot` recipe.
#'
#' @return The rgl device id, returned invisibly.
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
cell_plot_rgl <- function(object) {
  .validate_cell_plot(object)
  expect_rgl()
  expect_later()

  object$mapping$arrange <- NULL
  return(.render_cell_plot_rgl(build_cell_plot(object)))
}

#' Render a built cell plot with rgl
#'
#' Converts rendered colors, sizes, and alpha into rgl point properties.
#' Colors come from .cell_plot_rendered_colors() so builder-baked illumination
#' matches ggplot, Plotly, and base. Panel grids become an [rgl::layout3d()]
#' arrangement that keeps every row and column level, including empty panels.
#' Facet labels use dedicated gray strip regions; title and legend chrome use
#' reserved regions so placement does not depend on the upper-left panel.
#'
#' @param object A `cell_plot_built` object with a mapped `z` coordinate.
#'
#' @return The rgl device id, returned invisibly.
#'
#' @noRd
.render_cell_plot_rgl <- function(object) {
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

  device <- rgl::open3d(windowRect = c(100, 100, 1100, 1100))
  # Overlays are bitmaps tied to the viewport they were drawn in, so keep the
  # drawing calls around to repaint them whenever the window size changes.
  chrome <- list()

  title_id <- NULL
  legend_id <- NULL
  corner_id <- NULL
  col_strip_ids <- integer()
  row_strip_ids <- integer()
  resize_layout <- NULL
  if (!has_layout) {
    panel_ids <- rgl::subsceneInfo()$id
  } else {
    parent_id <- rgl::currentSubscene3d()
    layout_args <- list(
      n_row = n_row,
      n_col = n_col,
      need_col_strips = !is.null(facet_cols),
      need_row_strips = !is.null(facet_rows),
      need_title = need_title,
      need_legend = need_legend,
      need_subtitle = !is.null(subtitle)
    )
    layout <- do.call(
      .cell_rgl_facet_layout,
      c(
        layout_args,
        list(viewport = rgl::par3d("viewport", subscene = parent_id))
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
    if (!is.null(facet_cols) && !is.null(facet_rows)) {
      resize_layout <- function() {
        viewport <- rgl::par3d("viewport", subscene = parent_id)
        resized_layout <- do.call(
          .cell_rgl_facet_layout,
          c(layout_args, list(viewport = viewport))
        )
        .cell_rgl_set_layout_viewports(
          ids = ids,
          mat = resized_layout$mat,
          widths = resized_layout$widths,
          heights = resized_layout$heights,
          parent_viewport = viewport
        )
        return(invisible(NULL))
      }
    }
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
    .cell_rgl_set_listeners(
      c(panel_ids, col_strip_ids, corner_id, row_strip_ids, title_id, legend_id),
      panel_ids = panel_ids
    )
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
    chrome <- .cell_rgl_add_chrome(
      chrome,
      subscene = col_strip_ids[[col_index]],
      fn = .cell_rgl_strip_overlay,
      label = .cell_plot_facet_label(col_levels[[col_index]]),
      angle = 0,
      text_color = text_color,
      text_size = text_size,
      background_color = strip_background_color
    )
  }
  for (row_index in seq_along(row_strip_ids)) {
    chrome <- .cell_rgl_add_chrome(
      chrome,
      subscene = row_strip_ids[[row_index]],
      fn = .cell_rgl_strip_overlay,
      label = .cell_plot_facet_label(row_levels[[row_index]]),
      angle = .cell_plot_row_strip_angle,
      text_color = text_color,
      text_size = text_size,
      background_color = strip_background_color
    )
  }
  if (!is.null(corner_id)) {
    chrome <- .cell_rgl_add_chrome(
      chrome,
      subscene = corner_id,
      fn = .cell_rgl_strip_overlay,
      label = "",
      angle = 0,
      text_color = text_color,
      text_size = text_size,
      background_color = background_color
    )
  }

  if (!is.null(title_id)) {
    chrome <- .cell_rgl_add_chrome(
      chrome,
      subscene = title_id,
      fn = .cell_rgl_title_overlay,
      title = title,
      subtitle = subtitle,
      text_color = text_color,
      text_size = text_size,
      background_color = background_color
    )
  }
  if (!is.null(legend_id)) {
    chrome <- .cell_rgl_add_chrome(
      chrome,
      subscene = legend_id,
      fn = .cell_rgl_legend_overlay,
      text_color = text_color,
      text_size = text_size,
      background_color = background_color,
      categorical_legend = categorical_legend,
      continuous_legend = continuous_legend
    )
  }

  resources <- .cell_rgl_draw_chrome(chrome)
  .cell_rgl_register_chrome(
    device = device,
    chrome = chrome,
    resources = resources,
    resize = resize_layout
  )

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

#' Resize rgl layout subscenes
#'
#' Converts layout weights into pixel viewports and applies them to existing
#' subscenes. This keeps facet-strip dimensions synchronized when an rgl window
#' is resized without rebuilding the scenes or their contents.
#'
#' @param ids Subscene ids indexed by the integer ids in `mat`.
#' @param mat Integer layout matrix.
#' @param widths,heights Relative layout dimensions.
#' @param parent_viewport Four-value rgl parent viewport.
#'
#' @return `ids`, invisibly.
#'
#' @noRd
.cell_rgl_set_layout_viewports <- function(
  ids,
  mat,
  widths,
  heights,
  parent_viewport
) {
  parent_viewport <- as.numeric(parent_viewport)
  pixel_widths <- parent_viewport[[3L]] * widths / sum(widths)
  pixel_heights <- parent_viewport[[4L]] * heights / sum(heights)
  x_positions <- c(0, cumsum(pixel_widths))
  y_positions <- rev(c(0, cumsum(rev(pixel_heights))))[-1L]

  for (id_index in seq_along(ids)) {
    rows <- range(row(mat)[mat == id_index])
    cols <- range(col(mat)[mat == id_index])
    viewport <- c(
      x = x_positions[[cols[[1L]]]],
      y = y_positions[[rows[[2L]]]],
      width = sum(pixel_widths[seq.int(cols[[1L]], cols[[2L]])]),
      height = sum(pixel_heights[seq.int(rows[[1L]], rows[[2L]])])
    )
    rgl::par3d(viewport = viewport, subscene = ids[[id_index]])
  }
  return(invisible(ids))
}

#' Connect rgl subscene mouse actions to the data panels
#'
#' Data panels listen together, while chrome forwards its mouse actions to the
#' data panels without moving itself.
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

#' Record an overlay drawing call so it can be repainted later
#'
#' [rgl::bgplot3d()] rasterizes into a bitmap sized for the subscene viewport
#' it was drawn in, so overlays must be redrawn to survive a window resize.
#' Each entry stores the target subscene together with a zero-argument closure
#' that repaints it.
#'
#' @param chrome List of existing overlay entries.
#' @param subscene Subscene id the overlay belongs to.
#' @param fn Drawing function to call.
#' @param ... Arguments passed to `fn`.
#'
#' @return `chrome` with one entry appended.
#'
#' @noRd
.cell_rgl_add_chrome <- function(chrome, subscene, fn, ...) {
  force(fn)
  args <- list(...)
  entry <- list(
    subscene = subscene,
    draw = function() {
      return(invisible(do.call(fn, args)))
    }
  )
  return(c(chrome, list(entry)))
}

#' Paint every recorded overlay at the current window size
#'
#' Overlays are drawn into their own subscenes, so the current subscene is
#' restored afterwards to avoid leaking a chrome subscene to the caller.
#'
#' @param chrome List of overlay entries from .cell_rgl_add_chrome().
#'
#' @return A list with `object_ids` and `texture_files`, invisibly.
#'
#' @noRd
.cell_rgl_draw_chrome <- function(chrome) {
  resources <- list(object_ids = integer(), texture_files = character())
  if (length(chrome) == 0L) {
    return(invisible(resources))
  }
  current_subscene <- rgl::currentSubscene3d()
  finished <- FALSE
  on.exit(
    {
      rgl::useSubscene3d(current_subscene)
      if (!finished) {
        .cell_rgl_remove_chrome(resources)
      }
    },
    add = TRUE
  )
  for (entry in chrome) {
    rgl::useSubscene3d(entry$subscene)
    background_id <- as.integer(entry$draw())
    resources$object_ids <- unique(c(
      resources$object_ids,
      background_id,
      as.integer(rgl::rgl.attrib(background_id, "ids"))
    ))
    resources$texture_files <- unique(c(
      resources$texture_files,
      rgl::material3d("texture", id = background_id)
    ))
  }
  finished <- TRUE
  return(invisible(resources))
}

#' Remove resources belonging to old rgl overlays
#'
#' [rgl::bgplot3d()] creates both a background object and a textured quad.
#' Removing only the returned background id leaves the quad in the scene, so
#' all recorded object ids and their temporary texture files are discarded
#' together after a replacement overlay has been drawn.
#'
#' @param resources A list with `object_ids` and `texture_files`.
#'
#' @return `NULL`, invisibly.
#'
#' @noRd
.cell_rgl_remove_chrome <- function(resources) {
  if (length(resources$object_ids) > 0L) {
    live_ids <- as.integer(names(rgl::scene3d()$objects))
    object_ids <- intersect(as.integer(resources$object_ids), live_ids)
    if (length(object_ids) > 0L) {
      rgl::pop3d(id = object_ids)
    }
  }
  if (length(resources$texture_files) > 0L) {
    unlink(resources$texture_files)
  }
  return(invisible(NULL))
}

#' Drop one device from the rgl chrome registry
#'
#' Deletes its temporary texture files before removing the registry entry.
#' Scene objects need no cleanup here because this path is used when the
#' associated device has closed or has no overlays.
#'
#' @param key Character device id used as the registry key.
#'
#' @return `TRUE` when an entry was removed, `FALSE` otherwise, invisibly.
#'
#' @noRd
.cell_rgl_drop_chrome <- function(key) {
  if (!exists(key, envir = .cell_rgl_chrome_registry, inherits = FALSE)) {
    return(invisible(FALSE))
  }
  state <- get(key, envir = .cell_rgl_chrome_registry)
  unlink(state$resources$texture_files)
  rm(list = key, envir = .cell_rgl_chrome_registry)
  if (length(ls(envir = .cell_rgl_chrome_registry)) == 0L) {
    .cell_rgl_cancel_poll()
  }
  return(invisible(TRUE))
}

#' Registry of rgl devices with redrawable text overlays
#'
#' Keyed by device id, each element holds the overlay entries for that device
#' and the window rectangle they were last drawn at.
#'
#' @noRd
.cell_rgl_chrome_registry <- new.env(parent = emptyenv())

#' Seconds between idle checks of the rgl window size
#'
#' Dragging a window border produces no R command and no rgl event, so the
#' window size is polled while R is idle. The interval also acts as the
#' debounce window: a resize is repainted once the size holds still for one
#' full interval, so dragging repaints once at the end instead of continuously.
#'
#' @noRd
.cell_rgl_poll_interval <- 0.25

#' Whether an idle poll is already queued
#'
#' Guards against stacking [later::later()] callbacks, which would poll faster
#' and faster every time a plot is drawn.
#'
#' @noRd
.cell_rgl_poll_state <- new.env(parent = emptyenv())

#' Track a device so its overlays can be redrawn
#'
#' Registering with an empty `chrome` clears any previous registration.
#'
#' @param device rgl device id.
#' @param chrome List of overlay entries from .cell_rgl_add_chrome().
#' @param resources Current overlay object ids and texture files.
#' @param resize Optional zero-argument function that updates subscene geometry
#' after a window resize.
#'
#' @return `TRUE` if the device is now tracked, `FALSE` otherwise, invisibly.
#'
#' @noRd
.cell_rgl_register_chrome <- function(
  device,
  chrome,
  resources,
  resize = NULL
) {
  key <- as.character(device)
  .cell_rgl_drop_chrome(key)
  if (length(chrome) == 0L) {
    return(invisible(FALSE))
  }
  assign(
    key,
    list(
      chrome = chrome,
      resources = resources,
      resize = resize,
      window_rect = as.numeric(rgl::par3d("windowRect", dev = device))
    ),
    envir = .cell_rgl_chrome_registry
  )
  .cell_rgl_schedule_poll()
  return(invisible(TRUE))
}

#' Queue the next idle check of the rgl window size
#'
#' [later::later()] runs its callbacks while R is idle, including while the
#' console sits at the prompt, which is the only time a border drag can be
#' noticed. Polling stops once no device is tracked.
#'
#' @return `TRUE` if a poll was queued, `FALSE` otherwise, invisibly.
#'
#' @noRd
.cell_rgl_schedule_poll <- function() {
  if (isTRUE(.cell_rgl_poll_state$queued)) {
    return(invisible(FALSE))
  }
  if (length(ls(envir = .cell_rgl_chrome_registry)) == 0L) {
    .cell_rgl_cancel_poll()
    return(invisible(FALSE))
  }
  .cell_rgl_poll_state$queued <- TRUE
  .cell_rgl_poll_state$cancel <- later::later(
    function() {
      .cell_rgl_poll_state$queued <- FALSE
      .cell_rgl_poll_state$cancel <- NULL
      .cell_rgl_resize_poll()
      .cell_rgl_schedule_poll()
      return(invisible(NULL))
    },
    delay = .cell_rgl_poll_interval
  )
  return(invisible(TRUE))
}

#' Cancel the queued rgl resize poll
#'
#' [later::later()] returns a cancellation closure. Invoking it prevents a
#' queued idle callback from running after the last watched window closes or
#' the package unloads.
#'
#' @return `TRUE` when a callback was cancelled, `FALSE` otherwise, invisibly.
#'
#' @noRd
.cell_rgl_cancel_poll <- function() {
  cancel <- .cell_rgl_poll_state$cancel
  if (is.function(cancel)) {
    cancel()
  }
  .cell_rgl_poll_state$cancel <- NULL
  .cell_rgl_poll_state$queued <- FALSE
  return(invisible(is.function(cancel)))
}

#' Poll tracked devices and redraw overlays after settled resizes
#'
#' Closed devices are dropped. Changed devices are repainted only after their
#' window dimensions remain stable across two polls. An open device that fails
#' to repaint keeps its overlays and is retried on the next poll.
#'
#' @return `TRUE` while at least one device is still tracked.
#'
#' @noRd
.cell_rgl_resize_poll <- function() {
  keys <- ls(envir = .cell_rgl_chrome_registry)
  open_devices <- as.integer(rgl::rgl.dev.list())
  for (key in keys) {
    device <- as.integer(key)
    if (!device %in% open_devices) {
      .cell_rgl_drop_chrome(key)
      next
    }
    # A failed repaint must not stop resize handling for the other devices, and
    # must not discard the overlays an open window is still displaying, so the
    # device stays tracked and is retried on the next poll.
    tryCatch(
      .cell_rgl_check_device(key, device),
      error = function(condition) NULL
    )
  }
  keep_watching <- length(ls(envir = .cell_rgl_chrome_registry)) > 0L
  if (!keep_watching) {
    .cell_rgl_cancel_poll()
  }
  return(keep_watching)
}

#' Repaint one tracked device if its window size changed and has settled
#'
#' @param key Registry key for the device.
#' @param device rgl device id.
#'
#' @return `TRUE` if the device was repainted, `FALSE` otherwise, invisibly.
#'
#' @noRd
.cell_rgl_check_device <- function(key, device) {
  state <- get(key, envir = .cell_rgl_chrome_registry)
  window_rect <- as.numeric(rgl::par3d("windowRect", dev = device))
  if (isTRUE(all.equal(window_rect, state$window_rect))) {
    if (!is.null(state$pending_rect)) {
      state$pending_rect <- NULL
      assign(key, state, envir = .cell_rgl_chrome_registry)
    }
    return(invisible(FALSE))
  }
  if (!isTRUE(all.equal(window_rect, state$pending_rect))) {
    # Size is still changing, so the drag is not finished yet.
    state$pending_rect <- window_rect
    assign(key, state, envir = .cell_rgl_chrome_registry)
    return(invisible(FALSE))
  }
  .cell_rgl_redraw_chrome(device)
  return(invisible(TRUE))
}

#' Repaint the overlays of one tracked rgl device
#'
#' @param device rgl device id.
#'
#' @return `TRUE` if overlays were repainted, `FALSE` otherwise, invisibly.
#'
#' @noRd
.cell_rgl_redraw_chrome <- function(device) {
  key <- as.character(device)
  if (!exists(key, envir = .cell_rgl_chrome_registry, inherits = FALSE)) {
    return(invisible(FALSE))
  }
  if (!device %in% rgl::rgl.dev.list()) {
    .cell_rgl_drop_chrome(key)
    return(invisible(FALSE))
  }

  state <- get(key, envir = .cell_rgl_chrome_registry)
  previous_device <- as.integer(rgl::cur3d())
  if (previous_device != device) {
    rgl::set3d(device, silent = TRUE)
    on.exit(
      if (previous_device %in% rgl::rgl.dev.list()) {
        rgl::set3d(previous_device, silent = TRUE)
      },
      add = TRUE
    )
  }

  if (is.function(state$resize)) {
    state$resize()
  }
  resources <- .cell_rgl_draw_chrome(state$chrome)
  .cell_rgl_remove_chrome(state$resources)
  state$resources <- resources
  state$window_rect <- as.numeric(rgl::par3d("windowRect", dev = device))
  state$pending_rect <- NULL
  assign(key, state, envir = .cell_rgl_chrome_registry)
  return(invisible(TRUE))
}

#' Draw a facet label in a dedicated rgl strip
#'
#' The strip is a separate subscene from the data panel, so points cannot cover
#' its text. A gray background follows the default [ggplot2::facet_grid()]
#' appearance. Row strip labels are rotated vertically.
#'
#' @param label Facet level label.
#' @param angle Text rotation in degrees.
#' @param text_color,text_size,background_color Theme values.
#'
#' @return The value returned by [rgl::bgplot3d()], invisibly.
#'
#' @noRd
.cell_rgl_strip_overlay <- function(
  label,
  angle,
  text_color,
  text_size,
  background_color
) {
  return(invisible(rgl::bgplot3d(
    {
      graphics::par(
        mar = c(0, 0, 0, 0),
        bg = background_color,
        fg = text_color,
        cex = text_size / 11
      )
      graphics::plot(
        0,
        0,
        type = "n",
        xlim = c(0, 1),
        ylim = c(0, 1),
        axes = FALSE,
        xlab = "",
        ylab = "",
        xaxs = "i",
        yaxs = "i"
      )
      graphics::text(
        x = 0.5,
        y = 0.5,
        labels = label,
        srt = angle,
        col = text_color
      )
    },
    bg.color = background_color
  )))
}

#' Draw a title in a reserved rgl layout region
#'
#' @param title,subtitle Optional plot title and subtitle text.
#' @param text_color,text_size,background_color Theme values.
#'
#' @return The value returned by [rgl::bgplot3d()], invisibly.
#'
#' @noRd
.cell_rgl_title_overlay <- function(
  title,
  subtitle,
  text_color,
  text_size,
  background_color
) {
  return(invisible(rgl::bgplot3d(
    {
      graphics::par(
        mar = c(0, 1, 0, 1),
        bg = background_color,
        fg = text_color,
        cex = text_size / 11
      )
      graphics::plot(
        0,
        0,
        type = "n",
        xlim = c(0, 1),
        ylim = c(0, 1),
        axes = FALSE,
        xlab = "",
        ylab = "",
        xaxs = "i",
        yaxs = "i"
      )
      if (!is.null(title)) {
        graphics::text(
          x = 0.02,
          y = if (is.null(subtitle)) 0.5 else 0.65,
          labels = title,
          adj = c(0, 0.5),
          col = text_color,
          cex = 1.2,
          font = 2
        )
      }
      if (!is.null(subtitle)) {
        graphics::text(
          x = 0.02,
          y = if (is.null(title)) 0.5 else 0.25,
          labels = subtitle,
          adj = c(0, 0.5),
          col = text_color
        )
      }
    },
    bg.color = background_color
  )))
}

#' Draw a legend in a reserved rgl layout region
#'
#' @param text_color,text_size,background_color Theme values.
#' @param categorical_legend Optional list with `title`, `labels`, and `colors`.
#' @param continuous_legend Optional list with `title`, `limits`, and `colors`.
#'
#' @return The value returned by [rgl::bgplot3d()], invisibly.
#'
#' @noRd
.cell_rgl_legend_overlay <- function(
  text_color,
  text_size,
  background_color,
  categorical_legend = NULL,
  continuous_legend = NULL
) {
  return(invisible(rgl::bgplot3d(
    {
      graphics::par(
        mar = c(0, 0, 0, 0),
        bg = background_color,
        fg = text_color,
        col.axis = text_color,
        col.lab = text_color,
        cex = text_size / 11
      )
      if (!is.null(continuous_legend)) {
        .cell_rgl_draw_colorbar(continuous_legend, text_color = text_color)
      } else if (!is.null(categorical_legend)) {
        graphics::plot(
          0,
          0,
          type = "n",
          xlim = c(0, 1),
          ylim = c(0, 1),
          axes = FALSE,
          xlab = "",
          ylab = "",
          xaxs = "i",
          yaxs = "i"
        )
        graphics::legend(
          "center",
          legend = categorical_legend$labels,
          col = categorical_legend$colors,
          pch = 16,
          title = categorical_legend$title,
          bty = "n",
          text.col = text_color,
          title.col = text_color
        )
      }
    },
    bg.color = background_color
  )))
}

#' Draw a continuous colorbar for a numeric color scale
#'
#' Numeric mappings render as a continuous ramp rather than a discrete
#' `legend()` key, matching [cell_plot_interactive()]. The bar, ticks, and
#' title are positioned in relative coordinates so the bar keeps its
#' proportions in a narrow legend region.
#'
#' @param continuous_legend List with `title`, `limits`, and `colors`.
#' @param text_color Theme text color.
#'
#' @return `NULL`, invisibly.
#'
#' @noRd
.cell_rgl_draw_colorbar <- function(
  continuous_legend,
  text_color
) {
  n_stops <- max(100L, length(continuous_legend$colors))
  legend_colors <- grDevices::colorRampPalette(continuous_legend$colors)(
    n_stops
  )
  # image()/seq() need finite, strictly increasing limits. Constant or
  # all-non-finite color scales otherwise abort during faceted legend layout.
  limits <- continuous_legend$limits
  if (length(limits) < 2L || any(!is.finite(limits))) {
    limits <- c(0, 1)
  } else if (limits[1] == limits[2]) {
    pad <- max(abs(limits[1]) * 0.05, 1e-6)
    limits <- c(limits[1] - pad, limits[2] + pad)
  }
  continuous_legend$limits <- limits

  # A reserved legend region can be only a few dozen pixels wide. Base
  # graphics margins are measured in text lines, so they would consume the
  # whole region and collapse the bar into a thin strip with a clipped
  # title. Place the bar in relative coordinates instead.
  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot(
    0,
    0,
    type = "n",
    xlim = c(0, 1),
    ylim = c(0, 1),
    axes = FALSE,
    xlab = "",
    ylab = "",
    xaxs = "i",
    yaxs = "i"
  )
  bar_x <- c(0.18, 0.3)
  bar_y <- c(0.3, 0.7)
  title_y <- 0.77
  label_cex <- 0.85
  title_cex <- 0.9

  graphics::rasterImage(
    # Raster rows are painted from top to bottom, while scale colors are
    # ordered from low to high. Reverse them so high values appear at the top,
    # matching ggplot2's default vertical colorbar.
    grDevices::as.raster(matrix(rev(legend_colors), ncol = 1L)),
    xleft = bar_x[1],
    ybottom = bar_y[1],
    xright = bar_x[2],
    ytop = bar_y[2],
    interpolate = TRUE
  )

  tick_values <- pretty(limits, n = 4)
  tick_values <- tick_values[
    tick_values >= limits[1] & tick_values <= limits[2]
  ]
  tick_y <- stats::approx(x = limits, y = bar_y, xout = tick_values)$y
  tick_length <- diff(bar_x) * 0.2
  graphics::segments(
    bar_x[2],
    tick_y,
    bar_x[2] + tick_length,
    tick_y,
    col = text_color
  )
  graphics::text(
    x = bar_x[2] + tick_length * 1.4,
    y = tick_y,
    labels = tick_values,
    adj = c(0, 0.5),
    col = text_color,
    cex = label_cex
  )
  graphics::text(
    x = mean(bar_x),
    y = title_y,
    labels = continuous_legend$title,
    col = text_color,
    cex = title_cex
  )
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
