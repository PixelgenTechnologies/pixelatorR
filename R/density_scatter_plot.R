

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
    dens_grid$z[cbind(findInterval(x, dens_grid$x),
                      findInterval(y, dens_grid$y))]

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
#' @param marker1 Marker to plot along the x-axis.
#' @param marker2 Marker to plot along the y-axis.
#' @param facet_vars Variables to facet the plot by. If NULL, no faceting is done. If a character vector with 1 element,
#'                   the plot is faceted by rows. If a character vector with 2 elements, the plot is faceted by rows and
#'                   columns.
#' @param plot_gate A data.frame with columns 'xmin', 'xmax', 'ymin', 'ymax' to plot a gate. This data.frame can also
#'                  contain the variables in 'facet_vars' to plot different gates in different facets.
#' @param scale_density Scale the density to the maximum density per facet.
#' @param margin_density Add marginal density plots. Only supported when 'facet_vars' is NULL.
#' @param coord_fixed Fix the aspect ratio of the plot. Only supported when 'margin_density' is FALSE.
#' @param pt_size Size of the points.
#' @param alpha Transparency of the points.
#' @param layer Name of layer to plot. If NULL, the default layer is used.
#' @param grid_n Number of grid points to calculate the density.
#' @param colors Colors to use for the density plot. If NULL, the 'viridis' "turbo" palette is used.
#' @param ... Additional arguments to pass to 'MASS::kde2d'.
#'
#' @return A ggplot object.
#'
#' @import rlang
#'
#' @examples
#'
#' library(pixelatorR)
#' library(Seurat)
#'
#' # A mock-up Seurat Object
#' object <-
#'  CreateSeuratObject(counts = matrix(c(rpois(100000, 40),
#'                                       rpois(100000, 5))[sample(1:200000, 200000)],
#'                                     nrow = 100, ncol = 2000,
#'                                     dimnames = list(paste0("Feature", 1:100),
#'                                                     paste0("Cell", 1:2000))))
#'
#' object <-
#'   AddMetaData(object,
#'               metadata = data.frame(sample = rep(c("A", "B"), each = 1000),
#'                                     sample_type = rep(c("Unstimulated", "Stimulated"),
#'                                                       each = 500, times = 2),
#'                                     row.names = paste0("Cell", 1:2000)))
#'
#' plot_gate <-
#'   data.frame(xmin = c(70, 75),
#'              xmax = c(150, 155),
#'              ymin = c(50, 50),
#'              ymax = c(150, 150),
#'              sample = c("A", "B"))
#'
#'
#' DensityScatterPlot(object,
#'                    marker1 = "Feature1",
#'                    marker2 = "Feature2",
#'                    facet_vars = "sample",
#'                    plot_gate = plot_gate,
#'                    layer = "counts")
#'
#' @export
#'
DensityScatterPlot <- function(
    object,
    marker1,
    marker2,
    facet_vars = NULL,
    plot_gate = NULL,
    scale_density = TRUE,
    margin_density = FALSE,
    coord_fixed = TRUE,
    pt_size = 1,
    alpha = 1,
    layer = NULL,
    grid_n = 500,
    colors = NULL,
    ...
) {

  # Validate input parameters
  stopifnot(
    "'facet_vars' must either be NULL or be a character vector with 1 or 2 elements" =
      is.null(facet_vars) | length(facet_vars) %in% 1:2,
    "Variables in 'facet_vars' must be available in the objects 'meta.data' slot" =
      is.null(facet_vars) | all(facet_vars %in% colnames(object[[]])),

    "'marker1' must be a character of length 1" =
      is.character(marker1),
    "'marker2' must be a character of length 1" =
      is.character(marker2),
    "'marker1' must be available in the object" =
      marker1 %in% rownames(object),
    "'marker2' must be available in the object" =
      marker2 %in% rownames(object),

    "'grid_n' must be a numeric of length 1" =
      is.numeric(grid_n),
    "'scale_density' must be either TRUE or FALSE" =
      is.logical(scale_density),
    "'margin_density' must be either TRUE or FALSE" =
      is.logical(margin_density),

    "'pt_size' must be a numeric of length 1" =
      is.numeric(pt_size),
    "'alpha' must be a numeric of length 1" =
      is.numeric(alpha),
    "'layer' must be a character of length 1" =
      is.character(layer) | is.null(layer),
    "'coord_fixed' must be either TRUE or FALSE" =
      is.logical(coord_fixed),

    "'plot_gate' must be either NULL or a data.frame" =
      is.null(plot_gate) | is.data.frame(plot_gate),
    "'plot_gate' must have columns 'xmin', 'xmax', 'ymin', 'ymax'" =
      is.null(plot_gate) | all(c("xmin", "xmax", "ymin", "ymax") %in% colnames(plot_gate)),
    "'plot_gate' can't contain facetting variables that are not in 'facet_vars'" =
      is.null(plot_gate) | all(names(plot_gate) %in% c("xmin", "xmax", "ymin", "ymax", facet_vars)),

    "'object' must be a Seurat object" =
      inherits(object, "Seurat"),

    "Marginal density is not supported when 'facet_vars' is not NULL" =
      is.null(facet_vars) | !isTRUE(margin_density)
  )

  if(isTRUE(coord_fixed) && isTRUE(margin_density)) {
    warn("Fixed coordinates ('coord_fixed' = TRUE) is not supported when 'margin_density' is TRUE")
  }

  # Get data
  plot_data <-
    FetchData(object, vars = c(facet_vars, marker1, marker2),
              layer = layer)

  plot_data <-
    plot_data %>%
    rename(marker1 = !!marker1,
           marker2 = !!marker2) %>%
    group_by_at(facet_vars) %>%
    mutate(dens = .get2Ddensity(marker1, marker2, n = grid_n, ...)) %>%
    ungroup()

  if(isTRUE(scale_density)) {
    plot_data <-
      plot_data %>%
      group_by_at(facet_vars) %>%
      mutate(dens = dens / max(dens)) %>%
      ungroup()
  }

  # Set plot theme
  if(is.null(facet_vars)) {

    plot_theme <-
      theme_bw() +
      theme(panel.grid = element_blank(),
            panel.spacing = unit(0, "lines"),
            plot.margin = unit(c(0, 0, 0, 0), "lines"))

  } else {

    plot_theme <-
      theme_bw() +
      theme(panel.grid = element_blank())

  }

  plot_range <-
    range(c(plot_data$marker1, plot_data$marker2))

  if(!is.null(plot_gate)) {
    plot_range <-
      range(c(plot_data$marker1, plot_data$marker2,
              plot_gate$xmin, plot_gate$xmax,
              plot_gate$ymin, plot_gate$ymax))
  }

  # Make plot
  plot <-
    plot_data %>%
    ggplot(aes(marker1, marker2, color = dens)) +
    geom_hline(yintercept = 0, color = "gray") +
    geom_vline(xintercept = 0, color = "gray") +
    geom_point(size = pt_size,
               alpha = alpha,
               show.legend = FALSE) +
    scale_x_continuous(limits = plot_range) +
    scale_y_continuous(limits = plot_range) +
    plot_theme +
    labs(x = marker1,
         y = marker2)

  if(isTRUE(coord_fixed) && isFALSE(margin_density)) {
    plot <- plot + coord_fixed()
  }

  if(!is.null(colors)) {

    plot <-
      plot +
      scale_color_gradientn(colors = colors)

  } else {

    plot <-
      plot +
      scale_color_viridis_c(option = "H")
  }

  # Facet
  if(!is.null(facet_vars)) {

    if(length(facet_vars) == 1) {

      plot <-
        plot +
        facet_grid(rows = vars(!!!syms(facet_vars)))
    }

    if(length(facet_vars) == 2) {
      plot <-
        plot +
        facet_grid(rows = vars(!!!syms(facet_vars[1])),
                   cols = vars(!!!syms(facet_vars[2])))
    }

  }


  # Add gates
  if(!is.null(plot_gate)) {

    if(!is.null(facet_vars)) {

      join_vars <- intersect(facet_vars, colnames(plot_gate))

      if(length(join_vars) > 0) {
        gate_label <-
          plot_gate %>%
          left_join(plot_data,
                    by = join_vars)
      } else {
        gate_label <-
          plot_gate %>%
          cross_join(plot_data)
      }

    } else {
      gate_label <-
        plot_data %>%
        bind_cols(plot_gate)
    }

    gate_label <-
      gate_label %>%
      mutate(in_gate = marker1 > xmin & marker1 < xmax & marker2 > ymin & marker2 < ymax) %>%
      group_by_at(c(facet_vars, "in_gate", "xmin", "xmax", "ymin", "ymax")) %>%
      summarise(n = n(),
                .groups = "drop") %>%
      group_by_at(facet_vars) %>%
      mutate(p = 100 * n / sum(n)) %>%
      ungroup() %>%
      filter(in_gate) %>%
      mutate(label = paste(round(p, 1), "%"),
             x = (xmax + xmin) / 2,
             y = ymax)

    plot <-
      plot +
      geom_rect(data = plot_gate,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                inherit.aes = FALSE,
                fill = "transparent",
                color = "black",
                linetype = "dashed") +
      geom_text(data = gate_label,
                aes(x = x, y = y, label = label),
                inherit.aes = FALSE,
                color = "black",
                vjust = 1,
                size = 3)
  }

  # Add marginal density
  if(isTRUE(margin_density)) {

    density_plot_x <-
      plot_data %>%
      ggplot(aes(marker1)) +
      geom_density(fill = "lightgray",
                   color = NA) +
      scale_x_continuous(limits = plot_range) +
      scale_y_continuous(expand = expansion(0)) +
      theme_void() +
      theme(panel.spacing = unit(0, "lines"),
            plot.margin = unit(c(0, 0, 0, 0), "lines"))

    density_plot_y <-
      plot_data %>%
      ggplot(aes(marker2)) +
      geom_density(fill = "lightgray",
                   color = NA) +
      scale_x_continuous(limits = plot_range) +
      scale_y_continuous(expand = expansion(0)) +
      theme_void() +
      theme(panel.spacing = unit(0, "lines"),
            plot.margin = unit(c(0, 0, 0, 0), "lines")) +
      coord_flip()

    plot <-
      density_plot_x +
      plot_spacer() +
      plot +
      density_plot_y +
      plot_layout(widths = c(5, 1),
                  heights = c(1, 5))


  }


  return(plot)
}




