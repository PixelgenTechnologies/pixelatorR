#' @include generics.R
NULL

# Declarations used in package check
globalVariables(
  names = c('edges'),
  package = 'pixelatorR',
  add = TRUE
)

#' @param group_by A character specifying a column to group by
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
#' pxl_file <- system.file("extdata/five_cells",
#'                         "five_cells.pxl",
#'                         package = "pixelatorR")
#' seur_obj <- ReadMPX_Seurat(pxl_file)
#' seur_obj
#'
#' # Plot edge ranks with data.frame
#' EdgeRankPlot(seur_obj[[]])
#'
#' @export
#'
EdgeRankPlot.data.frame <- function (
  object,
  group_by = NULL,
  ...
) {

  # Check object
  stopifnot(
    "'object' must be a nonempty 'data.frame'-like object" =
      length(object) > 0,
    "column 'edges' is missing from 'object'" =
      "edges" %in% colnames(object)
  )
  stopifnot(
    "'edges' must be an integer vector" =
      inherits(object[["edges"]], what = "integer")
  )

  if (!is.null(group_by)) {
    stopifnot(
      "'group_by' must be a character of length 1" =
        is.character(group_by) &&
        (length(group_by) == 1)
    )
    if (!inherits(object[[group_by]], what = c("character", "factor"))) {
      abort(glue("Invalid class '{class(object[[group_by]])}' for column ",
                 "'{group_by}'. Expected a 'character' or 'factor'"))
    }
    if (!group_by %in% colnames(object)) {
      cli_alert_warning("'{group_by}' not found in 'object'")
    } else {
      object <- object %>% group_by_at(group_by)
    }
  }

  object <- object %>%
    mutate(rank = rank(-edges, ties.method = "random"))

  # Create edge rank plot
  edgerank_plot <-
    object %>%
      {
        if (!is.null(group_by)) {
          ggplot(., aes(rank, edges, color = !! sym(group_by)))
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
  group_by = NULL,
  ...
) {

  edgerank_plot <- EdgeRankPlot(object[[]], group_by = group_by)
  return(edgerank_plot)
}
