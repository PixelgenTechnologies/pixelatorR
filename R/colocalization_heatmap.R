# Declarations used in package check
globalVariables(
  names = c('estimate'),
  package = 'pixelatorR',
  add = TRUE
)

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
#' @param ... Parameters passed to \code{pheatmap}
#'
#' @return A heatmap
#'
#' @export
#'
ColocalizationHeatmap <- function (
    data,
    marker1_col = "marker_1",
    marker2_col = "marker_2",
    value_col = "estimate",
    colors = c("#053061", "#2166AC", "#4393C3", "#92C5DE",
               "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582",
               "#D6604D", "#B2182B", "#67001F"),
    return_plot_data = FALSE,
    symmetrise = TRUE,
    ...
) {

  # Check if pheatmap is installed
  expect_pheatmap()

  # Validate data
  if (!inherits(data, what = "tbl_df")) {
    abort(glue("'data' must be a 'tbl_df' object"))
  }
  stopifnot(
    "'marker_1' is not in 'data'" = marker1_col %in% names(data),
    "'marker_2' is not in 'data'" = marker2_col %in% names(data),
    "'estimate' is not in 'data'" = value_col %in% names(data),
    "'colors' must be a character vector with at least 2 colors" =
      inherits(colors, what = "character") &&
      (length(colors) > 1)
  )


  plot_data <-
    data %>%
    select(marker_1 = !!marker1_col,
           marker_2 = !!marker2_col,
           value = !!value_col)

  # Validate heatmap data
  if(nrow(distinct(select(plot_data, marker_1, marker_2))) != nrow(plot_data)) {
    abort("Invalid data format for a heatmap: There are multiple values for marker1 and marker2")
  }

  if(symmetrise) {
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

  if (return_plot_data) return(plot_data)

  # Set range for heatmap legend
  legend_range <- max(abs(plot_data))

  # Plot heatmap
  plot_data %>%
    pheatmap::pheatmap(
      breaks = seq(-legend_range, legend_range, length.out = 101),
      color = colorRampPalette(colors)(100),
      ...
    )
}
