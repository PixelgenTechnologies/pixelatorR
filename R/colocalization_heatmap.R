#' Plot DCA results
#'
#' Draws a heatmap of estimates obtained from a differential analysis
#' test computed with \code{\link{RunDCA}}.
#'
#' @param data A \code{tbl_df} object generated with \code{\link{RunDCA}}
#' @param marker1_col The name of the column with the first marker
#' @param marker2_col The name of the column with the second marker
#' @param value_col The name of the column with the value to plot, such as the estimate of a differential analysis test.
#' @param colors A character vector with colors
#' @param return_plot_data Return data formatted for plotting
#' instead of drawing the heatmap
#' @param symmetrise Set to \code{TRUE} if only the lower or upper triangle of marker combinations exist in \code{data},
#'                   each row will then be mirrored to fill the missing triangle of the heatmap.
#' @param legend_range A numeric vector of length 2 with the range of the legend. If NULL, the range is set to the
#'                     maximum absolute value of the data. If a value is outside this range, it is set to the closest
#'                     legend range limit.
#' @param legend_title The title of the legend
#' @param ... Parameters passed to \code{pheatmap}
#'
#' @return A heatmap
#'
#' @export
#'
ColocalizationHeatmap <- function(
  data,
  marker1_col = "marker_1",
  marker2_col = "marker_2",
  value_col = "estimate",
  colors = c(
    "#053061", "#2166AC", "#4393C3", "#92C5DE",
    "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582",
    "#D6604D", "#B2182B", "#67001F"
  ),
  return_plot_data = FALSE,
  symmetrise = TRUE,
  legend_range = NULL,
  legend_title = "",
  ...
) {
  # Check if pheatmap is installed
  expect_ComplexHeatmap()

  # Validate data
  if (!inherits(data, what = "tbl_df")) {
    abort(glue("'data' must be a 'tbl_df' object"))
  }
  stopifnot(
    "'marker1_col' is not in 'data'" =
      marker1_col %in% names(data),
    "'marker2_col' is not in 'data'" =
      marker2_col %in% names(data),
    "'value_col' is not in 'data'" =
      value_col %in% names(data),
    "'colors' must be a character vector with at least 2 colors" =
      inherits(colors, what = "character") &&
        (length(colors) > 1),
    "'return_plot_data' must be a logical" =
      inherits(return_plot_data, what = "logical"),
    "'symmetrise' must be a logical" =
      inherits(symmetrise, what = "logical"),
    "'legend_range' must be a numeric vector of length 2" =
      is.null(legend_range) || (inherits(legend_range, what = "numeric") && length(legend_range) == 2),
    "'legend_title' must be a character" =
      inherits(legend_title, what = "character") & length(legend_title) == 1
  )


  plot_data <-
    data %>%
    select(
      marker_1 = !!marker1_col,
      marker_2 = !!marker2_col,
      value = !!value_col
    )

  # Validate heatmap data
  if (nrow(distinct(select(plot_data, marker_1, marker_2))) != nrow(plot_data)) {
    abort("Invalid data format for a heatmap: There are multiple values for marker1 and marker2")
  }

  if (symmetrise) {
    # Symmetrise data
    plot_data <-
      plot_data %>%
      bind_rows(plot_data %>%
        rename(marker_1 = marker_2, marker_2 = marker_1)) %>%
      distinct()
  }

  # Pivot results
  plot_data <-
    plot_data %>%
    pivot_wider(names_from = "marker_2", values_from = "value", values_fill = 0) %>%
    column_to_rownames("marker_1")

  if (return_plot_data) {
    return(plot_data)
  }

  # Set range for heatmap legend
  if (is.null(legend_range)) {
    legend_range <- max(abs(plot_data)) * c(-1, 1)
  }

  plot_data[plot_data < legend_range[1]] <- legend_range[1]
  plot_data[plot_data > legend_range[2]] <- legend_range[2]

  # Plot heatmap
  plot_data %>%
    as.matrix() %>%
    ComplexHeatmap::pheatmap(
      breaks = seq(legend_range[1], legend_range[2], length.out = 101),
      color = colorRampPalette(colors)(100),
      heatmap_legend_param = list(title = legend_title),
      ...
    )
}
