#' Plot a colocalization heatmap
#'
#' Heatmap of a summary statistic between marker pairs (e.g. DCA from
#' \code{\link{RunDCA}}).
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
#' @param type Plot type: \code{"tiles"} for a filled-tile heatmap, or
#' \code{"dots"} for a dot plot where size maps to \code{size_col}. Dot plots
#' return a \code{ggplot} object.
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
#' @param highlight_pairs Optional \code{tbl_df}/data.frame with columns
#' \code{marker_1} and \code{marker_2} naming pairs to outline with a cell
#' border. Pairs whose markers are absent from \code{data} are ignored.
#' \code{NULL} disables highlighting. May also carry a column named by
#' \code{highlight_color_col} when mapping border colours for
#' \code{type = "dots"}.
#' @param highlight_colors Border colour(s) for highlighted cells. Single colour
#' when \code{highlight_color_col} is \code{NULL} (default \code{"black"};
#' tiles and dots). With \code{highlight_color_col} (\code{type = "dots"} only):
#' a named vector for discrete values, or two or more colours for numeric
#' values.
#' @param highlight_color_col Column in \code{highlight_pairs} used to colour
#' borders (\code{type = "dots"} only; errors with \code{"tiles"}).
#' \code{NULL} (default) draws a constant border colour.
#' @param highlight_stroke Border line width for highlighted cells.
#' @param highlight_shrink Fraction by which to reduce the width and height of
#' highlighted cell borders. Must be at least 0 and less than 1. The default
#' \code{0.1} leaves a small gap between adjacent highlights; \code{0} uses the
#' full cell size.
#' @param ... Parameters passed to \code{pheatmap}
#'
#' @return A \code{Heatmap} object/plot if \code{type = "tiles"} or a \code{ggplot}
#' object/plot if \code{type = "dots"}.
#'
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
#' # Outline selected marker pairs (constant colour)
#' ColocalizationHeatmap(
#'   dca_markers,
#'   type = "dots",
#'   size_range = c(1, 10),
#'   highlight_pairs = tibble(
#'     marker_1 = c("M1", "M3"),
#'     marker_2 = c("M2", "M5")
#'   ),
#'   highlight_colors = "black"
#' )
#'
#' # Map border colours from a column in highlight_pairs (dots only)
#' ColocalizationHeatmap(
#'   dca_markers,
#'   type = "dots",
#'   size_range = c(1, 10),
#'   highlight_pairs = tibble(
#'     marker_1 = c("M1", "M3"),
#'     marker_2 = c("M2", "M5"),
#'     database = c("string", "alphafold")
#'   ),
#'   highlight_color_col = "database",
#'   highlight_colors = c(string = "#C0392B", alphafold = "#6C3483"),
#'   highlight_shrink = 0.15
#' )
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
  highlight_pairs = NULL,
  highlight_colors = "black",
  highlight_color_col = NULL,
  highlight_stroke = 1.2,
  highlight_shrink = 0.1,
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
  assert_vector(highlight_colors, type = "character", n = 1)
  assert_single_value(highlight_color_col, "string", allow_null = TRUE)
  assert_single_value(highlight_stroke, type = "numeric")
  assert_single_value(highlight_shrink, type = "numeric")
  if (
    !is.finite(highlight_shrink) ||
      highlight_shrink < 0 ||
      highlight_shrink >= 1
  ) {
    cli::cli_abort(
      "{.arg highlight_shrink} must be at least 0 and less than 1."
    )
  }
  assert_class(size_range, "numeric")
  assert_length(size_range, 2)
  if (any(size_range[1] < 0)) {
    cli::cli_abort(
      c(
        "i" = "{.arg size_range} must have non-negative values",
        "x" = "Found {.val {sum(size_range < 0)}} negative values"
      )
    )
  }
  if (!is.null(legend_range)) {
    assert_class(legend_range, c("numeric", "integer"))
    assert_length(legend_range, 2)
  }
  assert_class(highlight_pairs, c("tbl_df", "data.frame"), allow_null = TRUE)
  if (!is.null(highlight_pairs)) {
    assert_col_in_data("marker_1", highlight_pairs)
    assert_col_in_data("marker_2", highlight_pairs)
  }
  if (is.null(highlight_color_col)) {
    if (length(highlight_colors) != 1) {
      cli::cli_abort(
        c(
          "i" = paste(
            "When {.arg highlight_color_col} is {.code NULL},",
            "{.arg highlight_colors} must be a single colour string."
          ),
          "x" = "You've supplied a vector of length {length(highlight_colors)}."
        )
      )
    }
  } else {
    if (type == "tiles") {
      cli::cli_abort(
        c(
          "x" = paste(
            "{.arg highlight_color_col} is only supported when",
            "{.arg type} = {.val dots}."
          ),
          "i" = paste(
            "Use {.code type = \"dots\"}, or omit {.arg highlight_color_col}",
            "for a constant {.arg highlight_colors}."
          )
        )
      )
    }
    if (is.null(highlight_pairs)) {
      cli::cli_abort(
        c(
          "x" = "{.arg highlight_color_col} requires {.arg highlight_pairs}.",
          "i" = "Pass pairs that include column {.val {highlight_color_col}}."
        )
      )
    }
    assert_col_in_data(highlight_color_col, highlight_pairs)
    highlight_col_vals <- highlight_pairs[[highlight_color_col]]
    if (is.numeric(highlight_col_vals)) {
      if (length(highlight_colors) < 2) {
        cli::cli_abort(
          c(
            "i" = paste(
              "When {.arg highlight_color_col} is numeric,",
              "{.arg highlight_colors} must have at least 2 colours",
              "for {.fn scale_color_gradientn}."
            ),
            "x" = "You've supplied a vector of length {length(highlight_colors)}."
          )
        )
      }
    } else if (is.character(highlight_col_vals) || is.factor(highlight_col_vals)) {
      levels_obs <- unique(as.character(highlight_col_vals))
      color_names <- names(highlight_colors)
      if (!is.null(color_names) && any(nzchar(color_names))) {
        missing_levels <- setdiff(levels_obs, color_names)
        if (length(missing_levels) > 0) {
          cli::cli_abort(
            c(
              "x" = paste(
                "{.arg highlight_colors} is missing named colours for",
                "{.arg highlight_color_col} level(s): {.val {missing_levels}}."
              ),
              "i" = "Provide a named colour for every observed level."
            )
          )
        }
      }
    } else {
      cli::cli_abort(
        c(
          "i" = paste(
            "{.arg highlight_color_col} must be numeric, character, or factor",
            "in {.arg highlight_pairs}."
          ),
          "x" = "Column {.val {highlight_color_col}} is {.cls {class(highlight_col_vals)}}."
        )
      )
    }
  }

  # Check if the data is grouped
  if (is.grouped_df(data)) {
    cli::cli_warn(
      c("i" = "The input data is grouped and will be ungrouped.")
    )
    data <- data %>% ungroup()
  }

  cols_keep <- c(marker1_col, marker2_col, value_col)
  if (type == "dots") {
    assert_col_in_data(size_col, data)
    cols_keep <- c(cols_keep, size_col)
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
      c(
        "i" = "Each row must represent a unique {.var {marker1_col}}/{.var {marker2_col}} pair",
        "x" = "Found duplicated pairs: {.val {dup_pairs}}"
      )
    )
  }

  if (symmetrise) {
    # Symmetrise data
    plot_data <-
      plot_data %>%
      bind_rows(plot_data %>%
        rename(
          !!sym(marker1_col) := !!sym(marker2_col),
          !!sym(marker2_col) := !!sym(marker1_col)
        )) %>%
      distinct()
  }

  # Keep only highlight pairs whose markers appear in the plot
  if (!is.null(highlight_pairs) && nrow(highlight_pairs) > 0) {
    plot_markers <- unique(c(
      as.character(plot_data[[marker1_col]]),
      as.character(plot_data[[marker2_col]])
    ))
    highlight_pairs <- highlight_pairs %>%
      mutate(
        marker_1 = as.character(marker_1),
        marker_2 = as.character(marker_2)
      ) %>%
      filter(
        marker_1 %in% plot_markers,
        marker_2 %in% plot_markers
      )
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
      pivot_wider(names_from = all_of(marker2_col), values_from = all_of(value_col), values_fill = 0) %>%
      column_to_rownames(marker1_col) %>%
      as.matrix()

    # Cluster rows for the dots method
    if (cluster_rows && (type == "dots")) {
      rows_clust <- dist(plot_data_wide, method = clustering_distance_rows) %>%
        hclust(method = clustering_method)
      plot_data <- plot_data %>%
        mutate(!!sym(marker1_col) := factor(!!sym(marker1_col), levels = with(rows_clust, labels[order])))
    }

    # Cluster columns for the dots method
    if (cluster_cols && (type == "dots")) {
      if (symmetrise && cluster_rows) {
        cols_clust <- rows_clust
      } else {
        cols_clust <- dist(t(plot_data_wide), method = clustering_distance_cols) %>%
          hclust(method = clustering_method)
      }
      plot_data <- plot_data %>%
        mutate(!!sym(marker2_col) := factor(!!sym(marker2_col), levels = with(cols_clust, labels[order])))
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
    )

    if (!is.null(highlight_pairs) && nrow(highlight_pairs) > 0) {
      highlight_keep <- c(marker1_col, marker2_col)
      if (!is.null(highlight_color_col)) {
        highlight_keep <- c(highlight_keep, highlight_color_col)
      }
      highlight_cells <- highlight_pairs %>%
        mutate(
          !!sym(marker1_col) := marker_1,
          !!sym(marker2_col) := marker_2
        ) %>%
        select(all_of(highlight_keep))
      if (symmetrise) {
        highlight_cells <- bind_rows(
          highlight_cells,
          highlight_cells %>%
            rename(
              !!sym(marker1_col) := !!sym(marker2_col),
              !!sym(marker2_col) := !!sym(marker1_col)
            )
        ) %>%
          distinct()
      }
      if (is.null(highlight_color_col)) {
        p <- p +
          geom_tile(
            data = highlight_cells,
            aes(!!sym(marker1_col), !!sym(marker2_col)),
            fill = NA,
            colour = highlight_colors,
            linewidth = highlight_stroke,
            width = 1 - highlight_shrink,
            height = 1 - highlight_shrink,
            inherit.aes = FALSE
          )
      } else {
        p <- p +
          geom_tile(
            data = highlight_cells,
            aes(
              !!sym(marker1_col),
              !!sym(marker2_col),
              colour = !!sym(highlight_color_col)
            ),
            fill = NA,
            linewidth = highlight_stroke,
            width = 1 - highlight_shrink,
            height = 1 - highlight_shrink,
            inherit.aes = FALSE
          )
        if (is.numeric(highlight_pairs[[highlight_color_col]])) {
          p <- p +
            scale_color_gradientn(
              colors = highlight_colors,
              name = highlight_color_col
            )
        } else {
          p <- p +
            scale_color_manual(
              values = highlight_colors,
              name = highlight_color_col
            )
        }
      }
    }

    p <- p +
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
    pheatmap_args <- list(...)
    if (!is.null(highlight_pairs) && nrow(highlight_pairs) > 0) {
      highlight_mat <- matrix(
        FALSE,
        nrow = nrow(plot_data_wide),
        ncol = ncol(plot_data_wide),
        dimnames = dimnames(plot_data_wide)
      )
      rn <- rownames(plot_data_wide)
      cn <- colnames(plot_data_wide)
      for (i in seq_len(nrow(highlight_pairs))) {
        a <- as.character(highlight_pairs$marker_1[i])
        b <- as.character(highlight_pairs$marker_2[i])
        if (a %in% rn && b %in% cn) {
          highlight_mat[a, b] <- TRUE
        }
        # Match dots: only mirror when symmetrise = TRUE
        if (symmetrise && b %in% rn && a %in% cn) {
          highlight_mat[b, a] <- TRUE
        }
      }

      user_cell_fun <- pheatmap_args$cell_fun
      pheatmap_args$cell_fun <- NULL
      cell_fun <- function(j, i, x, y, width, height, fill, ...) {
        if (!is.null(user_cell_fun)) {
          user_cell_fun(j, i, x, y, width, height, fill, ...)
        }
        if (isTRUE(highlight_mat[i, j])) {
          grid::grid.rect(
            x = x,
            y = y,
            width = width * (1 - highlight_shrink),
            height = height * (1 - highlight_shrink),
            gp = grid::gpar(
              col = highlight_colors,
              fill = NA,
              lwd = highlight_stroke
            )
          )
        }
      }

      p <- do.call(
        ComplexHeatmap::pheatmap,
        c(
          list(
            mat = plot_data_wide,
            breaks = seq(legend_range[1], legend_range[2], length.out = 101),
            color = colorRampPalette(colors)(100),
            heatmap_legend_param = list(title = legend_title),
            cluster_rows = cluster_rows,
            cluster_cols = cluster_cols,
            clustering_distance_rows = clustering_distance_rows,
            clustering_distance_cols = clustering_distance_cols,
            clustering_method = clustering_method,
            cell_fun = cell_fun
          ),
          pheatmap_args
        )
      )
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
  }

  return(p)
}
