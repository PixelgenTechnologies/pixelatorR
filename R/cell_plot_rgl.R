#' Render a cell plot with rgl
#'
#' Draws a [cell_plot()] recipe as an interactive 3D scatter with [rgl] and
#' returns it as an htmlwidget. The widget has no fixed size: it fills the
#' IDE viewer, a Quarto or R Markdown page, or any other container and
#' follows that container when it is resized.
#'
#' Drag to rotate and scroll to zoom. Panels in a grid share one camera and
#' move together. Grids support at most 10 rows and 20 columns.
#'
#' Titles, facet strips, and legends are HTML drawn over the WebGL canvas, so
#' they stay sharp at any resolution and their text can be selected. Their
#' size follows the theme text size in points. Numeric color mappings get a
#' colorbar and categorical mappings a discrete legend; legends show the
#' color scale without illumination.
#'
#' Node sizes are converted from relative units to point diameters in pixels,
#' and continuous sizes are grouped into at most 20 size bins. Occlusion
#' follows the camera. The `arrange` and `depth` mappings and
#' [cell_coord_rotate()] are ignored. Hover labels and interactive legend
#' filtering are not available.
#'
#' @param object A `cell_plot` recipe.
#'
#' @return An htmlwidget.
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

  object$mapping$arrange <- NULL
  built <- build_cell_plot(object)
  return(.cell_rgl_html_widget(built))
}

#' Draw a cell plot into an rgl htmlwidget
#'
#' Draws the built plot on a temporary null device, captures it with
#' [rgl::rglwidget()], and closes the device. The widget gets no width or
#' height so it fills its container; passing `NULL` explicitly also bypasses
#' the knitr figure-size default of [rgl::rglwidget()]. Chrome is added by
#' the render hook in `.cell_rgl_chrome_js`.
#'
#' @param object A `cell_plot_built` object.
#'
#' @return An rgl htmlwidget.
#'
#' @noRd
.cell_rgl_html_widget <- function(object) {
  canvas <- .cell_rgl_default_canvas_px
  previous <- as.integer(rgl::cur3d())
  device <- rgl::open3d(
    useNULL = TRUE,
    silent = TRUE,
    windowRect = c(0, 0, canvas, canvas)
  )
  on.exit(.cell_rgl_close_html_device(device, previous), add = TRUE)
  rendered <- .render_cell_plot_rgl(object, device = device)
  # snapshot = FALSE: a null device cannot draw the raster snapshot knitr
  # requests for PDF or Word output.
  widget <- rgl::rglwidget(
    width = NULL,
    height = NULL,
    snapshot = FALSE
  )
  if (!is.null(rendered$chrome)) {
    widget <- htmlwidgets::onRender(
      widget,
      htmlwidgets::JS(.cell_rgl_chrome_js),
      data = rendered$chrome
    )
  }
  return(widget)
}

#' JavaScript that draws chrome as HTML and fits the panels around it
#'
#' Runs after the widget renders and after every resize and canvas restart.
#' Adds a layer over the canvas with the title block, one strip per facet
#' level, and the legend, built from `.cell_rgl_chrome_spec()`. The browser
#' measures the title height, legend width, and strip thickness; the data
#' panel viewports fill the remaining rectangle. The layer ignores the
#' pointer while each element accepts it, so text is selectable and dragging
#' on chrome does not move the plot. Labels are set with `textContent`. A
#' canvas restart clears the widget element, so the layer is re-attached on
#' every layout. A later render removes the previous layer and reuses the
#' resize wrapper, so titles and strips are not stacked.
#'
#' @noRd
.cell_rgl_chrome_js <- "
function(el, x, data) {
  var rgl = el.rglinstance;
  if (!rgl) {
    return;
  }
  var theme = data.theme,
      block = 'position:absolute;box-sizing:border-box;pointer-events:auto;',
      make = function(parent, css, text) {
        var node = document.createElement('div');
        node.style.cssText = css;
        if (text !== null && text !== undefined) {
          node.textContent = String(text);
        }
        parent.appendChild(node);
        return node;
      },
      i;
  var previous = el.querySelector('.cell-rgl-chrome');
  if (previous) {
    previous.remove();
  }
  if (window.getComputedStyle(el).position === 'static') {
    el.style.position = 'relative';
  }
  var layer = make(el,
    'position:absolute;left:0;top:0;overflow:hidden;pointer-events:none;' +
    'font-family:sans-serif;line-height:1.2;color:' + theme.textColor + ';' +
    'font-size:' + theme.textSize + 'pt;');
  layer.className = 'cell-rgl-chrome';

  var title = null;
  if (data.title !== null || data.subtitle !== null) {
    title = make(layer, block + 'left:0;top:0;padding:0.5em 0.75em;white-space:nowrap;');
    if (data.title !== null) {
      make(title, 'font-size:1.2em;font-weight:bold;', data.title);
    }
    if (data.subtitle !== null) {
      make(title, '', data.subtitle);
    }
  }

  var stripCss = block + 'display:flex;align-items:center;justify-content:center;' +
        'overflow:hidden;font-size:0.8em;background:' + theme.stripBackgroundColor + ';',
      labelCss = 'flex-shrink:0;white-space:nowrap;',
      colStrips = [], rowStrips = [], colLabels = [], rowLabels = [], strip;
  for (i = 0; i < (data.colStrips || []).length; i++) {
    strip = make(layer, stripCss);
    colLabels.push(make(strip, labelCss + 'padding:0.5em;', data.colStrips[i]));
    colStrips.push(strip);
  }
  for (i = 0; i < (data.rowStrips || []).length; i++) {
    strip = make(layer, stripCss);
    rowLabels.push(make(strip,
      labelCss + 'padding:0.5em;writing-mode:vertical-rl;transform:rotate(180deg);',
      data.rowStrips[i]));
    rowStrips.push(strip);
  }

  var legend = null;
  if (data.legend !== null) {
    legend = make(layer, block + 'right:0;padding:0.75em;white-space:nowrap;');
    if (data.legend.title !== null && data.legend.title !== '') {
      make(legend, 'margin-bottom:0.5em;', data.legend.title);
    }
    if (data.legend.type === 'discrete') {
      for (i = 0; i < data.legend.labels.length; i++) {
        var row = make(legend, 'display:flex;align-items:center;font-size:0.8em;line-height:1.7;');
        make(row, 'flex-shrink:0;width:1em;height:1em;border-radius:50%;margin-right:0.5em;' +
          'background:' + data.legend.colors[i] + ';');
        make(row, '', data.legend.labels[i]);
      }
    } else {
      var colors = data.legend.colors,
          fill = colors.length > 1 ? 'linear-gradient(to top,' + colors.join(',') + ')' : colors[0],
          body = make(legend, 'display:flex;align-items:stretch;height:12em;padding:0.5em 0;'),
          ticks, tick, labels = [], labelWidth = 0;
      make(body, 'width:1.4em;background:' + fill + ';');
      ticks = make(body, 'position:relative;font-size:0.8em;');
      for (i = 0; i < data.legend.ticks.length; i++) {
        tick = data.legend.ticks[i];
        make(ticks, 'position:absolute;left:0;width:0.4em;height:1px;background:' +
          theme.textColor + ';bottom:' + (tick.at * 100) + '%;');
        labels.push(make(ticks,
          'position:absolute;left:0.7em;transform:translateY(50%);bottom:' + (tick.at * 100) + '%;',
          tick.label));
      }
      for (i = 0; i < labels.length; i++) {
        labelWidth = Math.max(labelWidth, labels[i].offsetWidth);
      }
      ticks.style.width = 'calc(0.7em + ' + labelWidth + 'px)';
    }
  }

  var layout = function() {
    var W = rgl.canvas.width,
        H = rgl.canvas.height,
        legendW = legend ? legend.offsetWidth : 0,
        titleH = 0, stripPx = 0, r, c, sub;
    if (layer.parentNode !== el) {
      el.appendChild(layer);
    }
    layer.style.width = W + 'px';
    layer.style.height = H + 'px';
    if (title) {
      title.style.width = Math.max(W - legendW, 0) + 'px';
      titleH = title.offsetHeight;
    }
    for (i = 0; i < colLabels.length; i++) {
      stripPx = Math.max(stripPx, colLabels[i].offsetHeight);
    }
    for (i = 0; i < rowLabels.length; i++) {
      stripPx = Math.max(stripPx, rowLabels[i].offsetWidth);
    }
    stripPx = Math.ceil(stripPx);
    var colStripH = colStrips.length ? stripPx : 0,
        rowStripW = rowStrips.length ? stripPx : 0,
        left = rowStripW,
        top = titleH + colStripH,
        panelW = Math.max(W - left - legendW, 1),
        panelH = Math.max(H - top, 1),
        cellW = panelW / data.nCol,
        cellH = panelH / data.nRow;
    for (i = 0; i < data.panels.length; i++) {
      sub = rgl.getObj(data.panels[i]);
      if (!sub || !sub.par3d || !sub.par3d.viewport) {
        continue;
      }
      r = Math.floor(i / data.nCol);
      c = i % data.nCol;
      sub.par3d.viewport.x = (left + c * cellW) / W;
      sub.par3d.viewport.y = (H - top - (r + 1) * cellH) / H;
      sub.par3d.viewport.width = cellW / W;
      sub.par3d.viewport.height = cellH / H;
    }
    for (c = 0; c < colStrips.length; c++) {
      colStrips[c].style.left = (left + c * cellW) + 'px';
      colStrips[c].style.top = titleH + 'px';
      colStrips[c].style.width = cellW + 'px';
      colStrips[c].style.height = colStripH + 'px';
    }
    for (r = 0; r < rowStrips.length; r++) {
      rowStrips[r].style.left = '0px';
      rowStrips[r].style.top = (top + r * cellH) + 'px';
      rowStrips[r].style.width = rowStripW + 'px';
      rowStrips[r].style.height = cellH + 'px';
    }
    if (legend) {
      legend.style.top = Math.max(top + (panelH - legend.offsetHeight) / 2, 0) + 'px';
    }
  };
  el.cellRglLayout = layout;
  if (!rgl.cellRglChrome) {
    var resize = rgl.resize,
        restart = rgl.restartCanvas;
    rgl.cellRglChrome = true;
    rgl.resize = function(element) {
      resize.call(this, element);
      if (el.cellRglLayout) {
        el.cellRglLayout();
      }
    };
    rgl.restartCanvas = function() {
      restart.call(this);
      if (el.cellRglLayout) {
        el.cellRglLayout();
      }
    };
  }
  layout();
  rgl.drawScene();
}
"

#' Close the null device used for an rgl htmlwidget
#'
#' Restores the previously current device when it is still open.
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

#' Size in pixels of the null device the scene is drawn on
#'
#' @noRd
.cell_rgl_default_canvas_px <- 1000L

#' Render a built cell plot with rgl
#'
#' Draws the data panels as an [rgl::layout3d()] grid that keeps every row
#' and column level, including empty panels, using the rendered colors from
#' .cell_plot_rendered_colors() so illumination matches the other renderers.
#' The root subscene carries the background color and ignores the pointer.
#' Titles, facet strips, and legends are returned as a `chrome` specification
#' for the HTML render hook.
#'
#' @param object A `cell_plot_built` object with a mapped `z` coordinate.
#' @param device An open rgl device to draw into.
#'
#' @return A list with `device`, the rgl device id, and `chrome`, the list
#' from `.cell_rgl_chrome_spec()` or `NULL` when the plot has no title,
#' strips, or legend.
#'
#' @noRd
.render_cell_plot_rgl <- function(object, device) {
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

  title <- object$annotation$title
  subtitle <- object$annotation$subtitle
  need_legend <- !is.null(mapping$color)
  has_chrome <- !is.null(facet_rows) ||
    !is.null(facet_cols) ||
    !is.null(title) ||
    !is.null(subtitle) ||
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

  rgl::set3d(device, silent = TRUE)

  if (!has_chrome) {
    panel_ids <- rgl::subsceneInfo()$id
  } else {
    parent_id <- rgl::currentSubscene3d()
    # With layout3d()'s default inherited mouse handling, setting the root
    # mouseMode to "none" below would also disable every data panel.
    panel_ids <- rgl::layout3d(
      matrix(seq_len(n_panels), nrow = n_row, ncol = n_col, byrow = TRUE),
      sharedMouse = FALSE,
      mouseMode = "replace"
    )
    panel_ids <- as.integer(panel_ids[seq_len(n_panels)])
    .cell_rgl_set_listeners(panel_ids, panel_ids = panel_ids)
    rgl::useSubscene3d(parent_id)
    rgl::bg3d(color = background_color)
    rgl::par3d(
      mouseMode = c(
        left = "none",
        right = "none",
        middle = "none",
        wheel = "none"
      ),
      subscene = parent_id
    )
  }

  for (panel_index in seq_len(n_panels)) {
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

  chrome <- NULL
  if (has_chrome) {
    chrome <- .cell_rgl_chrome_spec(
      panel_ids = panel_ids,
      n_row = n_row,
      n_col = n_col,
      title = title,
      subtitle = subtitle,
      col_labels = if (!is.null(facet_cols)) {
        vapply(col_levels, .cell_plot_facet_label, character(1))
      },
      row_labels = if (!is.null(facet_rows)) {
        vapply(row_levels, .cell_plot_facet_label, character(1))
      },
      categorical_legend = categorical_legend,
      continuous_legend = continuous_legend,
      theme = object$theme
    )
  }
  return(invisible(list(device = as.integer(device), chrome = chrome)))
}

#' Describe chrome for the HTML render hook
#'
#' Collects what `.cell_rgl_chrome_js` needs to draw the title block, facet
#' strips, and legend. Indexed vectors are wrapped in lists so they serialize
#' as JSON arrays even with one element, colors are converted to hexadecimal
#' because R color names are not CSS colors, and colorbar ticks are chosen
#' with [pretty()].
#'
#' @param panel_ids Data panel subscene ids in row-major order.
#' @param n_row,n_col Panel grid dimensions.
#' @param title,subtitle Plot title and subtitle, or `NULL`.
#' @param col_labels,row_labels Facet strip labels per column and row, or
#' `NULL` when that direction is not faceted.
#' @param categorical_legend Optional list with `title`, `labels`, and `colors`.
#' @param continuous_legend Optional list with `title`, `limits`, and `colors`.
#' @param theme The plot theme list.
#'
#' @return A list with `panels`, `nRow`, `nCol`, `title`, `subtitle`,
#' `colStrips`, `rowStrips`, `legend`, and `theme`. `legend` is `NULL` or a
#' list with `type` `"discrete"` (`title`, `labels`, `colors`) or
#' `"continuous"` (`title`, `colors`, `ticks`), where each tick has `at`, its
#' position along the bar from 0 at the bottom to 1 at the top, and `label`.
#'
#' @noRd
.cell_rgl_chrome_spec <- function(
  panel_ids,
  n_row,
  n_col,
  title,
  subtitle,
  col_labels,
  row_labels,
  categorical_legend,
  continuous_legend,
  theme
) {
  legend <- NULL
  if (!is.null(continuous_legend)) {
    colors <- continuous_legend$colors
    if (length(colors) == 0L) {
      colors <- "#000000"
    }
    legend <- list(
      type = "continuous",
      title = continuous_legend$title,
      colors = as.list(colors),
      ticks = .cell_rgl_colorbar_ticks(
        .cell_rgl_colorbar_limits(continuous_legend$limits)
      )
    )
  } else if (!is.null(categorical_legend)) {
    legend <- list(
      type = "discrete",
      title = categorical_legend$title,
      labels = as.list(as.character(categorical_legend$labels)),
      colors = as.list(categorical_legend$colors)
    )
  }
  return(list(
    panels = as.list(as.integer(panel_ids)),
    nRow = as.integer(n_row),
    nCol = as.integer(n_col),
    title = title,
    subtitle = subtitle,
    colStrips = if (!is.null(col_labels)) as.list(col_labels),
    rowStrips = if (!is.null(row_labels)) as.list(row_labels),
    legend = legend,
    theme = list(
      textSize = theme$text_size,
      textColor = .cell_colors_to_hex(theme$text_color),
      backgroundColor = .cell_colors_to_hex(theme$background_color),
      stripBackgroundColor = .cell_colors_to_hex(theme$strip_background_color)
    )
  ))
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

#' Tick marks for a continuous colorbar
#'
#' Uses [pretty()] breaks that fall inside the limits. Labels use the default
#' formatting, which switches to scientific notation for extreme values.
#'
#' @param limits Finite, strictly increasing length-two range.
#'
#' @return A list with one entry per tick holding `at`, the position along
#' the bar from 0 to 1, and `label`, the formatted value.
#'
#' @noRd
.cell_rgl_colorbar_ticks <- function(limits) {
  values <- pretty(limits, n = 4)
  values <- values[values >= limits[[1]] & values <= limits[[2]]]
  labels <- format(values, trim = TRUE)
  at <- (values - limits[[1]]) / diff(limits)
  return(lapply(
    seq_along(values),
    function(i) list(at = at[[i]], label = labels[[i]])
  ))
}

#' Connect rgl data panels so they share one camera
#'
#' Each data panel listens to every data panel.
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
#' Uses the same unlit `points` primitive as [rgl::plot3d()] `type = "p"`,
#' with `point_antialias` so WebGL draws round points. [rgl::points3d()]
#' accepts only a scalar point size, so nodes are split into a bounded number
#' of size groups while preserving per-node color and alpha vectors within
#' each group.
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
        alpha = alpha[idx],
        point_antialias = TRUE
      )
    )
  }
  return(invisible(ids))
}
