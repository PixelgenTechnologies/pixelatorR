#' @include generics.R
NULL

# Declarations used in package check
globalVariables(
  names = c('molecules'),
  package = 'pixelatorR',
  add = TRUE
)

#' @param group_by A character specifying a column to group by
#'
#' @rdname MoleculeRankPlot
#' @method MoleculeRankPlot data.frame
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
#' MoleculeRankPlot(seur_obj[[]])
#'
#' @export
#'
MoleculeRankPlot.data.frame <- function (
  object,
  group_by = NULL,
  ...
) {

  # Check object
  stopifnot(
    "'object' must be a nonempty 'data.frame'-like object" =
      length(object) > 0,
    "Either column 'edges' or 'molecules' must be present in 'object'" =
      any(c("edges", "molecules") %in% colnames(object))
  )

  molecules_column <-
    ifelse("molecules" %in% colnames(object), "molecules", "edges")

  if(!inherits(object[[molecules_column]], what = "integer")) glue("'{molecules_column}' must be an integer vector")

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

  object <-
    object %>%
    rename(molecules = any_of(molecules_column)) %>%
    mutate(rank = rank(-molecules, ties.method = "random"))

  # Create edge rank plot
  cellrank_plot <-
    object %>%
      {
        if (!is.null(group_by)) {
          ggplot(., aes(rank, molecules, color = !! sym(group_by)))
        } else {
          ggplot(., aes(rank, molecules))
        }
      } +
      geom_point(size = 0.5) +
      scale_x_log10() +
      scale_y_log10() +
      labs(x = "Component rank (by number of molecules)",
           y = "Number of molecules") +
      theme_minimal()

  return(cellrank_plot)
}

#' @rdname MoleculeRankPlot
#' @method MoleculeRankPlot Seurat
#' @concept plots
#'
#' @examples
#' library(pixelatorR)
#'
#' # Plot edge ranks with Seurat object
#' MoleculeRankPlot(seur_obj)
#'
#' @export
#'
MoleculeRankPlot.Seurat <- function (
  object,
  group_by = NULL,
  ...
) {

  cellrank_plot <- MoleculeRankPlot(object[[]], group_by = group_by)
  return(cellrank_plot)
}
