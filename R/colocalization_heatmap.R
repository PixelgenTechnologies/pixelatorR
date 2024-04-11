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
#' @param colors A character vector with colors
#' @param return_plot_data Return data formatted for plotting
#' instead of drawing the heatmap
#' @param ... Parameters passed to \code{pheatmap}
#'
#' @return A heatmap
#'
#' @export
#'
ColocalizationHeatmap <- function (
  data,
  colors = c("#053061", "#2166AC", "#4393C3", "#92C5DE",
             "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582",
             "#D6604D", "#B2182B", "#67001F"),
  return_plot_data = FALSE,
  ...
) {

  # Check if pheatmap is installed
  expect_pheatmap()

  # Validate data
  if (!inherits(data, what = "tbl_df")) {
    abort(glue("'data' must be a 'tbl_df' object"))
  }
  stopifnot(
    "'Invalid columns names'" = all((data %>% names()) ==
                                      c("estimate", "data_type", "target", "reference",
                                        "n1", "n2", "statistic", "p", "conf.low", "conf.high",
                                        "method", "alternative", "marker_1", "marker_2", "p_adj")),
    "'colors' must be a character vector with at least 2 colors" =
      inherits(colors, what = "character") &&
      (length(colors) > 1)
  )

  # Pivot results
  plot_data <-
    data %>%
    arrange(-estimate) %>%
    select(marker_1, marker_2, estimate) %>%
    bind_rows(select(., marker_2 = marker_1, marker_1 = marker_2, estimate)) %>%
    pivot_wider(names_from = "marker_2", values_from = "estimate", values_fill = 0) %>%
    column_to_rownames("marker_1")

  if (return_plot_data) return(plot_data)

  # Set range for heatmap legend
  legend_range <- max(abs(plot_data))

  # Plot heatmap
  plot_data %>%
    pheatmap::pheatmap(
      breaks = seq(-legend_range, legend_range, length.out = 100),
      color = colorRampPalette(colors)(100),
      ...
    )
}
