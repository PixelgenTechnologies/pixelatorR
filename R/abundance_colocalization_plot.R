#' Create an abundance/colocalization scatterplot
#'
#' Create a scatter plot of the abundance values for two sets of markers from a Seurat object. The points
#' (each corresponding to a cell) in the scatter plot are colored by the colocalization of the two markers.
#'
#' @param object A Seurat object.
#' @param markers_x,markers_y A character vector of markers to plot along the x- and y-axis respectively.
#' @param shared_scales Use the same scales for all plot panels.
#' @param coord_fixed Fix the aspect ratio of the plot.
#' @param pt_size Size or size range of the points.
#' @param draw_origo Draw lines through the origin (0, 0) in the plot.
#' @param layer Name of layer to fetch abundance data from. If NULL, the default layer is used.
#' @param coloc_score Name of the colocalization score to plot.
#' @param colors Colors to use for the colocalization score.
#'
#' @return A ggplot object.
#'
#' @examples
#'
#' library(pixelatorR)
#'
#' # Load example data as a Seurat object
#' pxl_file <- system.file(
#'   "extdata/five_cells",
#'   "five_cells.pxl",
#'   package = "pixelatorR"
#' )
#' seur_obj <- ReadMPX_Seurat(pxl_file)
#'
#' # Plot with data.frame
#' AbundanceColocalizationPlot(seur_obj, c("CD3E", "CD4"), c("CD19", "CD20"))
#'
#' @export
#'
AbundanceColocalizationPlot <- function(
  object,
  markers_x,
  markers_y,
  shared_scales = TRUE,
  coord_fixed = TRUE,
  pt_size = c(1, 4),
  draw_origo = TRUE,
  layer = NULL,
  coloc_score = "pearson_z",
  colors = c(
    "#1F395F", "#496389", "#728BB1", "#AABAD1", "gray90",
    "#E9AABF", "#CD6F8D", "#A23F5E", "#781534"
  )
) {
  # Validate input parameters
  stopifnot(
    "'object' must be a Seurat object" =
      inherits(object, "Seurat"),
    "'markers_x' must be a character vector of length > 0" =
      is.character(markers_x) & length(markers_x) > 0,
    "'markers_y' must be a character vector of length > 0" =
      is.character(markers_y) & length(markers_y) > 0,
    "All markers in 'markers_x' must be present in the Seurat object" =
      all(markers_x %in% rownames(object)),
    "All markers in 'markers_y' must be present in the Seurat object" =
      all(markers_y %in% rownames(object)),
    "'shared_scales' must be either TRUE or FALSE" =
      is.logical(shared_scales),
    "'coord_fixed' must be either TRUE or FALSE" =
      is.logical(coord_fixed),
    "'pt_size' must be a numeric of length 1 or 2" =
      is.numeric(pt_size) & length(pt_size) %in% c(1, 2),
    "'draw_origo' must be either TRUE or FALSE" =
      is.logical(draw_origo),
    "'colors' must be a vector of colors" =
      is.character(colors),
    "'coloc_score' must be a character of length 1" =
      is.character(coloc_score)
  )


  # Get data
  plot_data <-
    FetchData(object,
      vars = c(markers_x, markers_y),
      layer = layer
    ) %>%
    as_tibble(rownames = "component")

  plot_data <-
    plot_data %>%
    pivot_longer(
      cols = all_of(markers_x),
      names_to = "marker_x",
      values_to = "value_x"
    ) %>%
    pivot_longer(
      cols = all_of(markers_y),
      names_to = "marker_y",
      values_to = "value_y"
    )

  coloc_scores <-
    ColocalizationScores(object)

  if (!coloc_score %in% colnames(coloc_scores)) {
    abort(glue("'{coloc_score}' is not a valid colocalization score."))
  }

  if (!is.numeric(coloc_scores[[coloc_score]])) {
    abort(glue(
      "'{coloc_score}' must reference a column with numeric values, ",
      "got {class(coloc_scores[[coloc_score]])}."
    ))
  }

  coloc_scores <-
    coloc_scores %>%
    select(component, marker_1, marker_2, coloc_score = all_of(coloc_score)) %>%
    filter(marker_1 %in% c(markers_x, markers_y)) %>%
    filter(marker_2 %in% c(markers_x, markers_y)) %>%
    filter(marker_1 != marker_2) %>%
    bind_rows(select(., component, marker_1 = marker_2, marker_2 = marker_1, coloc_score))

  plot_data <-
    plot_data %>%
    left_join(coloc_scores,
      by = c("component",
        "marker_x" = "marker_1",
        "marker_y" = "marker_2"
      )
    ) %>%
    mutate(
      marker_x = factor(marker_x, levels = markers_x),
      marker_y = factor(marker_y, levels = markers_y)
    )

  # Make plot
  if (length(pt_size) == 1) {
    plot <-
      plot_data %>%
      ggplot(aes(value_x, value_y,
        color = coloc_score
      ))
  } else {
    plot <-
      plot_data %>%
      ggplot(aes(value_x, value_y,
        color = coloc_score,
        size = abs(coloc_score)
      )) +
      scale_size_continuous(range = pt_size)
  }

  if (draw_origo) {
    plot <- plot +
      geom_vline(xintercept = 0, color = "gray") +
      geom_hline(yintercept = 0, color = "gray")
  }

  plot <-
    plot +
    geom_point() +
    scale_color_gradientn(
      colors = colors,
      limits = max(abs(range(coloc_scores$coloc_score, na.rm = TRUE))) * c(-1, 1)
    ) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(
      x = "Abundance",
      y = "Abundance",
      color = "Colocalization Score",
      size = "|Colocalization Score|"
    )

  if (isTRUE(coord_fixed)) {
    plot <- plot + coord_fixed()
  }

  if (isTRUE(shared_scales)) {
    plot_range <-
      range(c(plot_data$value_x, plot_data$value_y))

    plot <-
      plot +
      facet_grid(marker_y ~ marker_x) +
      scale_x_continuous(limits = plot_range) +
      scale_y_continuous(limits = plot_range)
  } else {
    plot <-
      plot +
      facet_grid(marker_y ~ marker_x, scales = "free")
  }

  return(plot)
}
