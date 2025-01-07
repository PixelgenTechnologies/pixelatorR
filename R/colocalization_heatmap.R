#' Plot a colocalization heatmap
#'
#' Draws a heatmap of some summary statistic between marker pairs stored in a
#' \code{tbl_df}. A typical use case is to show the estimates of a differential
#' colocalization analysis test (\code{\link{RunDCA}}).
#'
#' @section Input:
#' The input data should be a \code{tbl_df} object with at least these three columns:
#'
#' 1. A column with the first marker name.
#' 2. A column with the second marker name.
#' 3. A column with the numeric value to plot, such as the estimate of a differential
#' analysis test or colocalization scores.
#'
#' Each marker pair can only appear once in the data. This means that if
#' you ran the test across multiple groups, you need to subset the data first.
#'
#' @param data A \code{tbl_df}
#' @param marker1_col The name of the column with the first marker
#' @param marker2_col The name of the column with the second marker
#' @param value_col The name of the column with numeric values to plot,
#' such as the estimate of a differential analysis test.
#' @param size_col The name of a column with numeric values to scale the
#' dots by. This is only used when \code{type = "dots"}.
#' @param size_col_transform A function to transform the values in \code{size_col}.
#' For instance, if the size_col values are p-values, you can use \code{\(x) {-log10(x)}}
#' to transform them to a more interpretative scale. Note that the size legend
#' label will be the value for \code{size_col} with a "_transformed" suffix.
#' The label can be changed manually after plotting. See examples below.
#' @param size_range A numeric vector of length 2 with dot size range.
#' @param colors A character vector with colors to create a color scale that
#' the \code{value_col} are mapped to.
#' @param cluster_rows,cluster_cols A logical value indicating whether to cluster
#' the rows and/or columns.
#' @param clustering_distance_rows,clustering_distance_cols The distance method
#' to use for clustering rows and columns. Can be any method accepted by \code{dist}.
#' @param clustering_method The clustering method to use. Can be any method
#' accepted by \code{hclust}.
#' @param type The type of plot to draw. Can be \code{"tiles"} for a typical
#' heatmap where each marker pair corresponds to a filled tile or \code{"dots"}
#' for a "dot plot". In the latter, the sizes of the dots are scaled by the
#' \code{size_col} column. This representation has the added advantage that
#' size differences can be used to highlight other important information,
#' such as significance. The dot plot is a \code{ggplot} object which can
#' be easily modified to customize the style.
#' @param return_plot_data Return data formatted for plotting instead of drawing
#' the heatmap.
#' @param symmetrise Set to \code{TRUE} if only the lower or upper triangle
#' of marker combinations exist in \code{data}, each row will then be
#' mirrored to fill the missing triangle of the heatmap.
#' @param legend_range A numeric vector of length 2 with the range of the
#' legend. If NULL, the range is set to the maximum absolute value of the
#' data. If a value is outside this range, it is set to the closest
#' legend range limit.
#' @param legend_title The title of the legend
#' @param ... Parameters passed to \code{pheatmap}
#'
#' @return A \code{Heatmap} object/plot if \code{type = "tiles"} or a \code{ggplot}
#' object/plot if \code{type = "dots"}.
#'
#' @examples
#' library(pixelatorR)
#' library(dplyr)
#' library(ggplot2)
#'
#' # Create a table with artificial DCA results
#' # for the example
#' set.seed(123)
#' dca_markers <- tidyr::expand_grid(
#'   marker_1 = paste0("M", 1:10),
#'   marker_2 = paste0("M", 1:10)
#' ) %>%
#'   filter(marker_1 > marker_2) %>%
#'   mutate(estimate = rnorm(n(), 0, 1), p_adj = 10^(-rnorm(n(), 3, 2)))
#' dca_markers
#'
#' # Typical heatmap
#' ColocalizationHeatmap(dca_markers)
#'
#' # Skip symmetrisation
#' ColocalizationHeatmap(dca_markers, symmetrise = FALSE)
#'
#' # Typical dot plot
#' # Note that the size legend title needs to be set manually
#' ColocalizationHeatmap(
#'   dca_markers,
#'   type = "dots",
#'   size_range = c(1, 10),
#'   size_col_transform = \(x) {
#'     -log10(x)
#'   }
#' ) &
#'   labs(size = "-log10p(p_adj)")
#'
#' # We can specify the order of the markers with factor levels
#' # but we need to turn off the clustering
#' ColocalizationHeatmap(
#'   dca_markers %>%
#'     mutate(
#'       marker_1 = factor(marker_1, levels = paste0("M", 1:10)),
#'       marker_2 = factor(marker_2, levels = paste0("M", 1:10))
#'     ),
#'   type = "dots",
#'   size_range = c(1, 10),
#'   size_col_transform = \(x) {
#'     -log10(x)
#'   },
#'   cluster_rows = FALSE,
#'   cluster_cols = FALSE
#' ) &
#'   labs(size = "-log10p(p_adj)")
#'
#' # If there are multiple comparisons per marker pair, we need to subset
#' # the data first
#' dca_markers_two_tests <- bind_rows(
#'   dca_markers %>% mutate(refrence = "ctrl", target = "stim1"),
#'   dca_markers %>% mutate(refrence = "ctrl", target = "stim2")
#' )
#' dca_markers_two_tests
#'
#' # This will now fail
#' \dontrun{
#' ColocalizationHeatmap(dca_markers_two_tests)
#' }
#'
#' # We need to subset the data first
#' ColocalizationHeatmap(
#'   dca_markers_two_tests %>%
#'     filter(target == "stim1")
#' )
#'
#' @export
#'
ColocalizationHeatmap <- function(
  data,
  marker1_col = "marker_1",
  marker2_col = "marker_2",
  value_col = "estimate",
  size_col = "p_adj",
  size_col_transform = NULL,
  size_range = c(0.1, 3),
  colors = c(
    "#053061", "#2166AC", "#4393C3", "#92C5DE",
    "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582",
    "#D6604D", "#B2182B", "#67001F"
  ),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  type = c("tiles", "dots"),
  return_plot_data = FALSE,
  symmetrise = TRUE,
  legend_range = NULL,
  legend_title = "",
  ...
) {
  # Check if pheatmap is installed
  expect_ComplexHeatmap()

  # Validate data
  assert_class(data, "tbl_df")
  type <- match.arg(type, choices = c("tiles", "dots"))
  assert_col_in_data(marker1_col, data)
  assert_col_in_data(marker2_col, data)
  assert_col_in_data(value_col, data)
  assert_function(size_col_transform, allow_null = TRUE)
  assert_vector(colors, type = "character")
  assert_single_value(return_plot_data, type = "bool")
  assert_single_value(symmetrise, type = "bool")
  assert_single_value(legend_title, "string")
  assert_class(size_range, "numeric")
  assert_length(size_range, 2)
  if (any(size_range[1] < 0)) {
    cli::cli_abort(
      c("i" = "{.arg size_range} must have non-negative values",
        "x" = "Found {.val {sum(size_range < 0)}} negative values")
    )
  }
  if (!is.null(legend_range)) {
    assert_class(legend_range, c("numeric", "integer"))
    assert_length(legend_range, 2)
  }

  cols_keep <- c(marker1_col, marker2_col, value_col)
  numeric_cols <- value_col
  if (type == "dots") {
    assert_col_in_data(size_col, data)
    cols_keep <- c(cols_keep, size_col)
    numeric_cols <- c(numeric_cols, size_col)
  }
  plot_data <-
    data %>%
    select(
      all_of(cols_keep)
    )

  # Validate data
  assert_col_class(marker1_col, plot_data, classes = c("character", "factor"))
  assert_col_class(marker2_col, plot_data, classes = c("character", "factor"))
  assert_col_class(value_col, plot_data, classes = "numeric")
  if (size_col %in% cols_keep) {
    assert_col_class(size_col, plot_data, classes = "numeric")
  }

  # Validate heatmap data
  dup_rows <- duplicated(select(plot_data, all_of(c(marker1_col, marker2_col))))
  if (sum(dup_rows) > 0) {
    dup_pairs <- select(plot_data, all_of(c(marker1_col, marker2_col))) %>%
      apply(1, paste, collapse = "/")
    cli::cli_abort(
      c("i" = "Each row must represent a unique {.var {marker1_col}}/{.var {marker2_col}} pair",
        "x" = "Found duplicated pairs: {.val {dup_pairs}}")
    )
  }

  if (symmetrise) {
    # Symmetrise data
    plot_data <-
      plot_data %>%
      bind_rows(plot_data %>%
        rename(marker_1 = marker_2, marker_2 = marker_1)) %>%
      distinct()
  }

  # Set range for heatmap legend
  if (is.null(legend_range)) {
    legend_range <- max(abs(plot_data %>% pull(all_of(value_col)))) * c(-1, 1)
  }

  # Cap values to legend range
  plot_data <- plot_data %>%
    mutate(across(all_of(value_col), ~ case_when(
      .x < legend_range[1] ~ legend_range[1],
      .x > legend_range[2] ~ legend_range[2],
      TRUE ~ .x
    )))

  # Transform size col
  if (!is.null(size_col_transform)) {
    plot_data <- plot_data %>%
      mutate(!!sym(size_col) := size_col_transform(!!sym(size_col)))
    size_col_label <- paste0(size_col, "_transformed")
  } else {
    size_col_label <- size_col
  }

  # Cast long table to a wide matrix
  # This is necessary for pheatmap (tiles) but not for the ggplot based version (dots)
  # However, we also do it conditionally for dots if clustering is enabled to get the
  # correct order of rows and columns
  if (type == "tiles" || (cluster_rows || cluster_cols)) {
    plot_data_wide <-
      plot_data %>%
      {
        if (!is.null(size_col) && (size_col %in% names(.))) {
          select(., -all_of(size_col))
        } else {
          .
        }
      } %>%
      pivot_wider(names_from = "marker_2", values_from = all_of(value_col), values_fill = 0) %>%
      column_to_rownames("marker_1") %>%
      as.matrix()

    # Cluster rows for the dots method
    if (cluster_rows && (type == "dots")) {
      rows_clust <- dist(plot_data_wide, method = clustering_distance_rows) %>%
        hclust(method = clustering_method)
      plot_data <- plot_data %>%
        mutate(marker_1 = factor(marker_1, levels = with(rows_clust, labels[order])))
    }

    # Cluster columns for the dots method
    if (cluster_cols && (type == "dots")) {
      cols_clust <- dist(t(plot_data_wide), method = clustering_distance_cols) %>%
        hclust(method = clustering_method)
      plot_data <- plot_data %>%
        mutate(marker_2 = factor(marker_2, levels = with(cols_clust, labels[order])))
    }
  }

  # Return plot data if requested
  if (return_plot_data) {
    if (type == "dots") {
      return(plot_data)
    } else {
      return(plot_data_wide)
    }
  }

  # Plot heatmap
  if (type == "dots") {
    p <- ggplot(
      plot_data,
      aes(!!sym(marker1_col),
        !!sym(marker2_col),
        fill = !!sym(value_col),
        size = !!sym(size_col)
      )
    ) +
      geom_point(shape = 21) +
      scale_size(range = size_range) +
      scale_y_discrete(limits = rev) +
      scale_x_discrete(position = "top") +
      scale_fill_gradientn(colors = colors, limits = legend_range) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 0)) +
      coord_fixed() +
      labs(size = size_col_label)
  } else {
    p <- plot_data_wide %>%
      ComplexHeatmap::pheatmap(
        breaks = seq(legend_range[1], legend_range[2], length.out = 101),
        color = colorRampPalette(colors)(100),
        heatmap_legend_param = list(title = legend_title),
        cluster_rows = cluster_rows,
        cluster_cols = cluster_cols,
        clustering_distance_rows = clustering_distance_rows,
        clustering_distance_cols = clustering_distance_cols,
        clustering_method = clustering_method,
        ...
      )
  }

  return(p)
}
