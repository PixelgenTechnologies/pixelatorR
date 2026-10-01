#' Render a built cell plot with ggplot2
#'
#' Creates the default static 2D representation of a cell plot. Three-
#' dimensional layouts use the mapped x and y coordinates. Depth is mapped
#' to `z` by default and is represented through draw order and apparent point
#' size. Resolved hex colors, including
#' build-stage illumination, are drawn directly, while the unilluminated color
#' scale metadata supplies the legend.
#'
#' @param object A `cell_plot_built` object.
#' @param limits Optional shared x and y ranges used by animation frames.
#'
#' @return A `ggplot` object.
#'
#' @importFrom rlang .data
#'
#' @noRd
.render_cell_plot_ggplot <- function(object, limits = NULL) {
  pixelatorR:::assert_class(object, "cell_plot_built", arg = "object")

  plot_data <- object$data
  relative_size <- .cell_plot_projected_sizes(object)
  plot_data$.cell_size <- .cell_relative_size_to_ggplot(relative_size)
  plot_data$.cell_alpha <- object$alpha$resolved
  plot_data$.cell_color <- I(.cell_plot_rendered_colors(object))

  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data[[object$mapping$x]],
      y = .data[[object$mapping$y]],
      color = .data$.cell_color,
      size = .data$.cell_size,
      alpha = .data$.cell_alpha
    )
  ) +
    ggplot2::geom_point()

  if (!is.null(object$mapping$color)) {
    # AsIs colors bypass scale transformation while the scale metadata below
    # remains available to construct the legend.
    if (object$color$type == "continuous") {
      plot <- plot +
        ggplot2::scale_color_gradientn(
          colours = object$color$colors,
          limits = object$color$limits,
          oob = scales::squish,
          na.value = object$color$na_color,
          name = .cell_plot_legend_title(object)
        )
    } else {
      plot <- plot +
        ggplot2::scale_color_manual(
          values = object$color$colors,
          limits = object$color$limits,
          na.value = object$color$na_color,
          name = .cell_plot_legend_title(object),
          drop = FALSE
        )
    }
  }

  plot <- plot +
    ggplot2::scale_size_identity() +
    ggplot2::scale_alpha_identity()
  if (is.null(limits)) {
    plot <- plot + ggplot2::coord_fixed()
  } else {
    plot <- plot +
      ggplot2::coord_fixed(
        xlim = limits$x,
        ylim = limits$y,
        expand = FALSE
      )
  }
  plot <- plot +
    ggplot2::theme_void(base_size = object$theme$text_size) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(
        fill = object$theme$background_color,
        color = NA
      ),
      panel.background = ggplot2::element_rect(
        fill = object$theme$background_color,
        color = NA
      ),
      text = ggplot2::element_text(color = object$theme$text_color),
      plot.title = ggplot2::element_text(color = object$theme$text_color),
      plot.subtitle = ggplot2::element_text(color = object$theme$text_color),
      strip.text = ggplot2::element_text(
        color = object$theme$text_color,
        margin = .cell_plot_facet_strip_margin
      ),
      strip.text.y.left = ggplot2::element_text(
        angle = .cell_plot_row_strip_angle
      ),
      strip.background = ggplot2::element_rect(
        fill = object$theme$strip_background_color,
        color = NA
      ),
      legend.text = ggplot2::element_text(color = object$theme$text_color),
      legend.title = ggplot2::element_text(color = object$theme$text_color)
    ) +
    ggplot2::labs(
      title = object$annotation$title,
      subtitle = object$annotation$subtitle
    )

  if (!is.null(object$grid)) {
    facet_rows <- if (is.null(object$grid$rows)) {
      ggplot2::vars()
    } else {
      ggplot2::vars(!!rlang::sym(object$grid$rows))
    }
    facet_cols <- if (is.null(object$grid$cols)) {
      ggplot2::vars()
    } else {
      ggplot2::vars(!!rlang::sym(object$grid$cols))
    }
    plot <- plot +
      ggplot2::facet_grid(
        rows = facet_rows,
        cols = facet_cols,
        drop = FALSE,
        switch = "y"
      )
  }

  return(plot)
}
