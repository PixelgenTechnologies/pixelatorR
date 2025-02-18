#' Get density of points in 2 dimensions.
#'
#' Calculate the kernel density of points in 2 dimensions.
#'
#' @param x A numeric vector.
#' @param y A numeric vector.
#' @param n Create a square n by n grid to compute density.
#' @param ... Additional arguments to pass to 'MASS::kde2d'.
#' @return The density within each square.
#'
#' @noRd
#'
#' @examples
#'
#' library(pixelatorR)
#'
#' x <- rnorm(1000)
#' y <- rnorm(1000)
#'
#' get_density(x, y, n = 100)
#'
.get2Ddensity <- function(
  x,
  y,
  n = 500,
  ...
) {
  # Check that MASS is installed
  expect_MASS()

  dens_grid <-
    MASS::kde2d(x, y, n = n, ...)

  dens <-
    dens_grid$z[
      cbind(
        findInterval(x, dens_grid$x),
        findInterval(y, dens_grid$y)
      )
    ]

  return(dens)
}

#' Create a density scatter / pseudocolor plot.
#'
#' Create a density scatter plot, also known as a pseudocolor plot, of two markers from a Seurat object. The points in
#' the scatter plot are colored by the point density, determined by 2D kernel density estimation. Optionally, the plot
#' can be faceted by one or two variables and a gate can be plotted.
#'
#'
#' @param object A Seurat object.
#' @param marker1 Name of first marker to plot.
#' @param marker2 Name of second marker to plot.
#' @param facet_vars Optional character vector of length 1 or 2 specifying variables to facet by.
#' @param plot_gate Optional data frame containing gate parameters. If provided, `gate_type` must also be specified.
#' @param gate_type Optional string specifying the gate type. Must be either "rectangle" or "quadrant".
#' Required if `plot_gate` is provided.
#' @param grid_n Number of grid points in each direction to compute density.
#' @param scale_density Whether to scale the density within each facet.
#' @param margin_density Whether to plot density plots in the margins.
#' @param pt_size Point size.
#' @param alpha Point transparency.
#' @param layer Optional string specifying a layer to plot.
#' @param coord_fixed Whether to use fixed coordinate ratio.
#' @param annotation_params Optional list of parameters to pass to `geom_text` for gate annotations. Common parameters
#' include color (text color), vjust (vertical justification), hjust (horizontal justification), and size (text size).
#' @param colors Optional character vector of colors to use for the color scale.
#' @param ... Additional arguments to pass to `MASS::kde2d`.
#'
#' @return A ggplot object.
#'
#'
#' @examples
#'
#' library(pixelatorR)
#' library(Seurat)
#'
#' set.seed(123)
#'
#' # A mock-up Seurat Object
#' object <-
#'   CreateSeuratObject(counts = matrix(
#'     c(
#'       rpois(100000, 40),
#'       rpois(100000, 5)
#'     )[sample(1:200000, 200000)],
#'     nrow = 100, ncol = 2000,
#'     dimnames = list(
#'       paste0("Feature", 1:100),
#'       paste0("Cell", 1:2000)
#'     )
#'   ))
#'
#' object <-
#'   AddMetaData(object,
#'     metadata = data.frame(
#'       sample = rep(c("A", "B"), each = 1000),
#'       sample_type = rep(c("Unstimulated", "Stimulated"),
#'         each = 500, times = 2
#'       ),
#'       row.names = paste0("Cell", 1:2000)
#'     )
#'   )
#'
#' plot_gate <-
#'   data.frame(
#'     xmin = c(20, 20),
#'     xmax = c(70, 70),
#'     ymin = c(20, 20),
#'     ymax = c(60, 60),
#'     sample = c("A", "B")
#'   )
#'
#'
#' DensityScatterPlot(object,
#'   marker1 = "Feature1",
#'   marker2 = "Feature2",
#'   facet_vars = "sample",
#'   plot_gate = plot_gate,
#'   gate_type = "rectangle",
#'   layer = "counts"
#' )
#'
#'
#' # Create a Seurat object with random data
#' counts <- matrix(rpois(200, 5), nrow = 10)
#' rownames(counts) <- paste0("Feature", 1:10)
#' colnames(counts) <- paste0("Cell", 1:20)
#' object <- CreateSeuratObject(counts = counts)
#'
#' # Create a basic density scatter plot
#' DensityScatterPlot(
#'   object,
#'   marker1 = "Feature1",
#'   marker2 = "Feature2"
#' )
#'
#' # Create a density scatter plot with a quadrant gate
#' # Define quadrant gate parameters
#' quad_gate <- data.frame(
#'   x = 2, # x-coordinate for vertical line
#'   y = 3 # y-coordinate for horizontal line
#' )
#'
#' # Plot with quadrant gate and custom annotations
#' DensityScatterPlot(
#'   object,
#'   marker1 = "Feature1",
#'   marker2 = "Feature2",
#'   plot_gate = quad_gate,
#'   gate_type = "quadrant",
#'   annotation_params = list(
#'     color = "darkblue",
#'     size = 4,
#'     fontface = "bold"
#'   )
#' )
#'
#' @export
#'
DensityScatterPlot <- function(
  object,
  marker1,
  marker2,
  facet_vars = NULL,
  plot_gate = NULL,
  gate_type = c("rectangle", "quadrant"),
  grid_n = 500,
  scale_density = TRUE,
  margin_density = TRUE,
  pt_size = 1,
  alpha = 1,
  layer = NULL,
  coord_fixed = TRUE,
  annotation_params = NULL,
  colors = NULL,
  ...
) {
  # Validate input parameters
  assert_class(object, "Seurat")
  assert_vector(facet_vars, type = "character", n = 1, allow_null = TRUE)
  assert_max_length(facet_vars, n = 2, allow_null = TRUE)
  for (facet_var in facet_vars) {
    assert_col_in_data(facet_var, object[[]], allow_null = TRUE)
  }
  assert_single_value(marker1, type = "string")
  assert_single_value(marker2, type = "string")
  assert_x_in_y(marker1, rownames(object))
  assert_x_in_y(marker2, rownames(object))
  assert_single_value(grid_n, type = "integer")
  assert_single_value(scale_density, type = "bool")
  assert_single_value(margin_density, type = "bool")
  assert_single_value(pt_size, type = "numeric")
  assert_single_value(alpha, type = "numeric")
  assert_single_value(layer, type = "string", allow_null = TRUE)
  assert_single_value(coord_fixed, type = "bool")
  if (!is.null(annotation_params)) {
    assert_class(annotation_params, "list")
  }
  if (!is.null(plot_gate)) {
    assert_class(plot_gate, "data.frame")
    if (is.null(gate_type)) {
      cli::cli_abort(
        "{.code gate_type} must be specified when {.code plot_gate} is provided"
      )
    }
    gate_type <- match.arg(gate_type)
    if (gate_type == "rectangle") {
      assert_x_in_y(c("xmin", "xmax", "ymin", "ymax"), colnames(plot_gate))
    } else {
      assert_x_in_y(c("x", "y"), colnames(plot_gate))
    }
    allowed_cols <- if (gate_type == "rectangle") {
      c("xmin", "xmax", "ymin", "ymax")
    } else {
      c("x", "y")
    }
    check <- !names(plot_gate) %in% c(allowed_cols, facet_vars)
    if (sum(check) > 0) {
      cli::cli_abort(
        c(
          "i" = "{.code plot_gate} can only contain gate-specific columns and faceting variables",
          "x" = "The following columns are not allowed: ",
          " " = "{.val {names(plot_gate)[check]}}"
        )
      )
    }
  }

  # Set margin_density to FALSE if facet_vars is provided
  if (!is.null(facet_vars)) {
    if (isTRUE(margin_density)) {
      cli::cli_alert_warning(
        "Marginal density ({.codemargin_density = TRUE}) is not supported with faceting. Setting {.code margin_density = FALSE}"
      )
      margin_density <- FALSE
    }
  }

  if (isTRUE(coord_fixed) && isTRUE(margin_density)) {
    coord_fixed <- FALSE
    cli::cli_alert_warning(
      "Setting {.code coord_fixed = FALSE} as it is not compatible with {.code margin_density = TRUE}"
    )
  }

  plot_data <-
    FetchData(object, vars = c(facet_vars, marker1, marker2), layer = layer)

  plot_data <-
    plot_data %>%
      rename(
        marker1 = !!marker1,
        marker2 = !!marker2
      ) %>%
      group_by_at(if (is.null(facet_vars)) character(0) else facet_vars) %>%
      mutate(dens = .get2Ddensity(marker1, marker2, n = grid_n, ...)) %>%
      ungroup()

  if (isTRUE(scale_density)) {
    plot_data <-
      plot_data %>%
        group_by_at(if (is.null(facet_vars)) character(0) else facet_vars) %>%
        mutate(dens = dens / max(dens)) %>%
        ungroup()
  }

  # Set plot theme
  if (is.null(facet_vars)) {
    plot_theme <-
      theme_bw() +
        theme(
          panel.grid = element_blank(),
          panel.spacing = unit(0, "lines"),
          plot.margin = unit(c(0, 0, 0, 0), "lines")
        )
  } else {
    plot_theme <-
      theme_bw() +
        theme(panel.grid = element_blank())
  }

  plot_range <-
    range(c(plot_data$marker1, plot_data$marker2))

  if (!is.null(plot_gate)) {
    if (gate_type == "rectangle") {
      plot_range <- range(
        c(
          plot_data$marker1,
          plot_data$marker2,
          plot_gate$xmin,
          plot_gate$xmax,
          plot_gate$ymin,
          plot_gate$ymax
        )
      )
    } else if (gate_type == "quadrant") {
      plot_range <- range(
        c(
          plot_data$marker1,
          plot_data$marker2,
          plot_gate$x,
          plot_gate$y
        )
      )
    }
  }

  # Make plot
  plot <-
    plot_data %>%
      ggplot(aes(marker1, marker2, color = dens)) +
      geom_hline(yintercept = 0, color = "gray") +
      geom_vline(xintercept = 0, color = "gray") +
      geom_point(
        size = pt_size,
        alpha = alpha,
        show.legend = FALSE
      ) +
      scale_x_continuous(limits = plot_range) +
      scale_y_continuous(limits = plot_range) +
      plot_theme +
      labs(
        x = marker1,
        y = marker2
      )

  if (isTRUE(coord_fixed) && isFALSE(margin_density)) {
    plot <- plot + coord_fixed()
  }

  # Add colors
  if (!is.null(colors)) {
    plot <-
      plot +
        scale_color_gradientn(colors = colors)
  } else {
    plot <-
      plot +
        scale_color_viridis_c(option = "turbo")
  }

  # Facet
  if (!is.null(facet_vars)) {
    if (length(facet_vars) == 1) {
      plot <-
        plot +
          facet_grid(rows = vars(!!!syms(facet_vars)))
    }

    if (length(facet_vars) == 2) {
      plot <-
        plot +
          facet_grid(
            rows = vars(!!!syms(facet_vars[1])),
            cols = vars(!!!syms(facet_vars[2]))
          )
    }
  }

  # Add gates
  if (!is.null(plot_gate)) {
    if (gate_type == "rectangle") {
      if (!is.null(facet_vars)) {
        join_vars <- intersect(facet_vars, colnames(plot_gate))
        if (length(join_vars) > 0) {
          gate_label <- plot_gate %>% left_join(plot_data, by = join_vars)
        } else {
          gate_label <- plot_gate %>% cross_join(plot_data)
        }
      } else {
        gate_label <- plot_data %>% cross_join(plot_gate)
      }

      gate_label <- gate_label %>%
        mutate(
          in_gate = marker1 > xmin &
            marker1 < xmax &
            marker2 > ymin &
            marker2 < ymax
        ) %>%
        group_by_at(
          c(facet_vars, "in_gate", "xmin", "xmax", "ymin", "ymax")
        ) %>%
        summarise(n = n(), .groups = "drop") %>%
        group_by_at(c(facet_vars, "xmin", "xmax", "ymin", "ymax")) %>%
        mutate(p = 100 * n / sum(n)) %>%
        ungroup() %>%
        filter(in_gate) %>%
        mutate(
          label = paste(round(p, 1), "%"),
          x = (xmin + xmax) / 2,
          y = ymax
        )

      plot <- plot +
        geom_rect(
          data = plot_gate,
          aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
          inherit.aes = FALSE,
          fill = "transparent",
          color = "black",
          linetype = "dashed"
        ) +
        do.call(
          geom_text,
          c(
            list(
              data = gate_label,
              mapping = aes(x = x, y = y, label = label),
              inherit.aes = FALSE
            ),
            if (is.null(annotation_params)) {
              list(color = "black", vjust = 1, size = 3)
            } else {
              annotation_params
            }
          )
        )
    } else if (gate_type == "quadrant") {
      assert_x_in_y(c("x", "y"), colnames(plot_gate))

      if (!is.null(facet_vars)) {
        join_vars <- intersect(facet_vars, colnames(plot_gate))
        if (length(join_vars) > 0) {
          gate_label <- plot_gate %>% left_join(plot_data, by = join_vars)
        } else {
          gate_label <- plot_gate %>% cross_join(plot_data)
        }
      } else {
        gate_label <- plot_data %>% cross_join(plot_gate)
      }

      x_plot_range <- range(plot_data$marker1)
      y_plot_range <- range(plot_data$marker2)
      y_offset <- diff(y_plot_range) * 0.05

      gate_label <- gate_label %>%
        mutate(
          quadrant = case_when(
            marker1 < x & marker2 > y ~ "top_left",
            marker1 >= x & marker2 > y ~ "top_right",
            marker1 < x & marker2 <= y ~ "bottom_left",
            marker1 >= x & marker2 <= y ~ "bottom_right"
          )
        ) %>%
        group_by_at(c(facet_vars, "x", "y", "quadrant")) %>%
        summarise(n = n(), .groups = "drop") %>%
        group_by_at(c(facet_vars, "x", "y")) %>%
        mutate(p = 100 * n / sum(n)) %>%
        ungroup() %>%
        mutate(
          label = paste(round(p, 1), "%"),
          # Position labels in inner corners
          x_label = case_when(
            quadrant == "bottom_left" ~ x_plot_range[1],
            quadrant == "top_left" ~ x_plot_range[1],
            quadrant == "bottom_right" ~ x_plot_range[2],
            quadrant == "top_right" ~ x_plot_range[2]
          ),
          y_label = case_when(
            quadrant == "bottom_left" ~ y - y_offset,
            quadrant == "top_left" ~ y_plot_range[2],
            quadrant == "bottom_right" ~ y - y_offset,
            quadrant == "top_right" ~ y_plot_range[2]
          ),
          # Set text alignment based on quadrant
          hjust = case_when(
            quadrant %in% c("top_left", "bottom_left") ~ 0,
            quadrant %in% c("top_right", "bottom_right") ~ 1
          ),
          vjust = case_when(
            quadrant %in% c("top_left", "top_right") ~ 0,
            quadrant %in% c("bottom_left", "bottom_right") ~ 1
          )
        )

      # Add quadrant lines
      plot <- plot +
        geom_vline(
          xintercept = unique(plot_gate$x),
          linetype = "dashed",
          color = "black"
        ) +
        geom_hline(
          yintercept = unique(plot_gate$y),
          linetype = "dashed",
          color = "black"
        ) +
        do.call(
          geom_text,
          c(
            list(
              data = gate_label,
              mapping = aes(
                x = x_label,
                y = y_label,
                label = label,
                hjust = hjust,
                vjust = vjust
              ),
              inherit.aes = FALSE,
              check_overlap = TRUE
            ),
            if (is.null(annotation_params)) {
              list(color = "black", size = 4)
            } else {
              annotation_params
            }
          )
        )
    } else {
      cli::cli_abort("Unknown gate type provided to plot_gate")
    }
  }

  # Add marginal density
  if (isTRUE(margin_density)) {
    density_plot_x <-
      plot_data %>%
        ggplot(aes(marker1)) +
        geom_density(
          fill = "lightgray",
          color = NA
        ) +
        scale_x_continuous(limits = plot_range) +
        scale_y_continuous(expand = expansion(0)) +
        theme_void() +
        theme(
          panel.spacing = unit(0, "lines"),
          plot.margin = unit(c(0, 0, 0, 0), "lines")
        )

    density_plot_y <-
      plot_data %>%
        ggplot(aes(marker2)) +
        geom_density(
          fill = "lightgray",
          color = NA
        ) +
        scale_x_continuous(limits = plot_range) +
        scale_y_continuous(expand = expansion(0)) +
        theme_void() +
        theme(
          panel.spacing = unit(0, "lines"),
          plot.margin = unit(c(0, 0, 0, 0), "lines")
        ) +
        coord_flip()

    plot <-
      density_plot_x +
        plot_spacer() +
        plot +
        density_plot_y +
        plot_layout(
          widths = c(5, 1),
          heights = c(1, 5)
        )
  }

  return(plot)
}
