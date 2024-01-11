#' @include generics.R
NULL


#' @param group.by A character specifying a column to group by
#'
#' @import dplyr
#' @import ggplot2
#' @import cli
#' @import glue
#'
#' @rdname EdgeRankPlot
#' @method EdgeRankPlot data.frame
#' @concept plots
#'
#' @examples
#'
#' library(pixelatorR)
#'
#' # Load example data as a Seurat object
#' pxl_file <- system.file("extdata/PBMC_10_cells",
#'                         "Sample01_test.pxl",
#'                         package = "pixelatorR")
#' seur_obj <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)
#' seur_obj
#'
#' # Plot edge ranks with data.frame
#' EdgeRankPlot(seur_obj[[]])
#'
#' @export
#'
EdgeRankPlot.data.frame <- function (
  object,
  group.by = NULL,
  ...
) {

  # Check object
  stopifnot("'object' must be a nonempty 'data.frame'-like object" = length(object) > 0,
            "column 'edges' is missing from 'object'" = "edges" %in% colnames(object))
  stopifnot("'edges' must be an integer vector" = inherits(object[["edges"]], what = "integer"))

  if (!is.null(group.by)) {
    stopifnot(
      "'group.by' must be a character of length 1" =
        is.character(group.by) &&
        (length(group.by) == 1)
    )
    if (!inherits(object[[group.by]], what = c("character", "factor"))) {
      abort(glue("Invalid class '{class(object[[group.by]])}' for column ",
                 "'{group.by}'. Expected a 'character' or 'factor'"))
    }
    if (!group.by %in% colnames(object)) {
      cli_alert_warning("'{group.by}' not found in 'object'")
    } else {
      object <- object %>% group_by_at(group.by)
    }
  }

  object <- object %>%
    mutate(rank = rank(-edges, ties.method = "random"))

  # Create edge rank plot
  edgerank_plot <-
    object %>%
      {
        if (!is.null(group.by)) {
          ggplot(., aes(rank, edges, color = !! sym(group.by)))
        } else {
          ggplot(., aes(rank, edges))
        }
      } +
      geom_point(size = 0.5) +
      scale_x_log10() +
      scale_y_log10() +
      labs(x = "Component rank (by number of edges)",
           y = "Number of edges") +
      theme_minimal()

  return(edgerank_plot)
}

#' @rdname EdgeRankPlot
#' @method EdgeRankPlot Seurat
#' @concept plots
#'
#' @examples
#' library(pixelatorR)
#'
#' # Plot edge ranks with Seurat object
#' EdgeRankPlot(seur_obj)
#'
#' @export
#'
EdgeRankPlot.Seurat <- function (
  object,
  group.by = NULL,
  ...
) {

  edgerank_plot <- EdgeRankPlot(object[[]], group.by = group.by)
  return(edgerank_plot)
}
