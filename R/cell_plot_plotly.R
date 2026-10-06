#' Render a built cell plot with Plotly
#'
#' The Plotly backend of [cell_plot_interactive()]. Converts resolved colors,
#' sizes, and alpha into Plotly marker properties. Panel grids become a
#' Cartesian arrangement of scenes that keeps every row and column level,
#' including empty panels.
#'
#' @param object A `cell_plot_built` object with a mapped `z` coordinate.
#'
#' @return A Plotly htmlwidget.
#'
#' @noRd
.render_cell_plot_plotly <- function(object) {
  assert_class(object, "cell_plot_built", arg = "object")

  mapping <- object$mapping
  categorical <- identical(object$color$type, "categorical")
  n_nodes <- nrow(object$data)

  plot_data <- object$data
  # Plotly 3D markers take one opacity per trace, so node alpha travels
  # inside the marker color instead.
  plot_data$.cell_color <- plotly::toRGB(
    rep_len(.cell_plot_rendered_colors(object), n_nodes),
    alpha = rep_len(object$alpha$resolved, n_nodes)
  )
  # Depth sizing is ggplot-only; Plotly uses resolved sizes as-is.
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
        "x" = "{.fn cell_plot_interactive} supports at most 10 facet rows.",
        "i" = "{.arg rows} creates {length(row_levels)} facet rows."
      )
    )
  }
  if (length(col_levels) > 20) {
    cli::cli_abort(
      c(
        "x" = "{.fn cell_plot_interactive} supports at most 20 facet columns.",
        "i" = "{.arg cols} creates {length(col_levels)} facet columns."
      )
    )
  }
  panel_gap <- if (is.null(object$grid)) c(0, 0) else c(0.02, 0.04)

  plot <- NULL
  layout <- list()
  annotations <- list()
  panel_index <- 0L

  for (row_index in seq_along(row_levels)) {
    for (col_index in seq_along(col_levels)) {
      panel_index <- panel_index + 1L
      scene <- paste0("scene", if (panel_index == 1L) "" else panel_index)
      domain <- list(
        x = c(
          (col_index - 1) / length(col_levels) + panel_gap[1],
          col_index / length(col_levels) - panel_gap[1]
        ),
        y = c(
          1 - row_index / length(row_levels) + panel_gap[2],
          1 - (row_index - 1) / length(row_levels) - panel_gap[2]
        )
      )

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

      if (nrow(panel_data) == 0) {
        # Empty panels still need a trace so their scene is drawn, matching
        # facet_grid(drop = FALSE).
        plot <- .cell_plotly_trace(
          plot,
          scene = scene,
          x = mean(ranges$x),
          y = mean(ranges$y),
          z = mean(ranges$z),
          marker = list(size = 0, opacity = 0, line = list(width = 0)),
          hoverinfo = "skip",
          showlegend = FALSE
        )
      } else {
        # One visible trace per panel, with mixed categorical colors and
        # NA / out-of-limit points. Nodes keep their input order because
        # scatter3d resolves occlusion against the scene camera.
        plot <- .cell_plotly_trace(
          plot,
          scene = scene,
          x = panel_data[[mapping$x]],
          y = panel_data[[mapping$y]],
          z = panel_data[[mapping$z]],
          marker = list(
            color = panel_data$.cell_color,
            size = panel_data$.cell_size,
            line = list(width = 0)
          ),
          showlegend = FALSE
        )
      }

      # Legend entries come from legend-only traces so that every level is
      # listed, including levels without nodes. Dummy coordinates are required
      # because Plotly.js omits empty traces from the legend.
      if (panel_index == 1L && categorical) {
        for (level in object$color$limits) {
          plot <- .cell_plotly_trace(
            plot,
            scene = scene,
            x = ranges$x[1],
            y = ranges$y[1],
            z = ranges$z[1],
            marker = list(
              color = .cell_colors_to_hex(unname(object$color$colors[level])),
              size = 8,
              line = list(width = 0)
            ),
            name = level,
            legendgroup = level,
            showlegend = TRUE,
            visible = "legendonly",
            hoverinfo = "skip"
          )
        }
      }

      if (panel_index == 1L && !categorical && !is.null(mapping$color)) {
        palette <- .cell_colors_to_hex(unname(object$color$colors))
        stops <- if (length(palette) == 1) {
          list(list(0, palette), list(1, palette))
        } else {
          lapply(seq_along(palette), function(stop_index) {
            list((stop_index - 1) / (length(palette) - 1), palette[stop_index])
          })
        }
        plot <- .cell_plotly_trace(
          plot,
          scene = scene,
          x = ranges$x[1],
          y = ranges$y[1],
          z = ranges$z[1],
          marker = list(
            color = object$color$limits[1],
            cmin = object$color$limits[1],
            cmax = object$color$limits[2],
            colorscale = stops,
            showscale = TRUE,
            colorbar = list(
              title = list(text = .cell_plot_legend_title(object)),
              outlinewidth = 0
            ),
            size = 0,
            opacity = 0,
            line = list(width = 0)
          ),
          hoverinfo = "skip",
          showlegend = FALSE
        )
      }

      layout[[scene]] <- list(
        domain = domain,
        xaxis = list(range = ranges$x, visible = FALSE),
        yaxis = list(range = ranges$y, visible = FALSE),
        zaxis = list(range = ranges$z, visible = FALSE),
        aspectmode = "data",
        bgcolor = object$theme$background_color
      )

      if (!is.null(object$grid)) {
        label <- paste(
          c(
            if (!is.null(facet_rows)) row_levels[[row_index]],
            if (!is.null(facet_cols)) col_levels[[col_index]]
          ),
          collapse = " | "
        )
        annotations <- c(annotations, list(list(
          text = label,
          x = mean(domain$x),
          y = domain$y[2],
          xref = "paper",
          yref = "paper",
          xanchor = "center",
          yanchor = "bottom",
          showarrow = FALSE,
          bgcolor = object$theme$strip_background_color,
          font = list(
            color = object$theme$text_color,
            size = object$theme$text_size
          )
        )))
      }
    }
  }

  title <- object$annotation$title
  if (!is.null(object$annotation$subtitle)) {
    title <- paste(
      c(title, paste0("<sup>", object$annotation$subtitle, "</sup>")),
      collapse = "<br>"
    )
  }

  layout$paper_bgcolor <- object$theme$background_color
  layout$plot_bgcolor <- object$theme$background_color
  layout$font <- list(
    color = object$theme$text_color,
    size = object$theme$text_size
  )
  layout$title <- if (!is.null(title)) {
    list(text = title, x = 0, xanchor = "left")
  }
  layout$annotations <- annotations
  layout$margin <- list(
    t = if (is.null(title)) 40 else 80,
    r = 80,
    b = 40,
    l = 40
  )
  layout$showlegend <- categorical
  layout$legend <- list(
    title = list(text = .cell_plot_legend_title(object) %||% "")
  )
  layout$hovermode <- "closest"

  return(do.call(plotly::layout, c(list(p = plot), layout)))
}

#' Add one marker trace to a Plotly cell plot
#'
#' Attaches a 3D scatter trace to the scene of one panel.
#'
#' @param plot A Plotly object.
#' @param scene The scene identifier of the target panel.
#' @param x,y,z Marker coordinates.
#' @param ... Further trace arguments passed to [plotly::add_trace()].
#'
#' @return The updated Plotly object.
#'
#' @noRd
.cell_plotly_trace <- function(plot, scene, x, y, z, ...) {
  trace <- rlang::list2(
    type = "scatter3d",
    mode = "markers",
    x = x,
    y = y,
    z = z,
    scene = scene,
    ...
  )
  trace <- Filter(Negate(is.null), trace)

  if (is.null(plot)) {
    return(do.call(plotly::plot_ly, trace))
  }

  return(do.call(
    plotly::add_trace,
    c(list(p = plot, inherit = FALSE), trace)
  ))
}
