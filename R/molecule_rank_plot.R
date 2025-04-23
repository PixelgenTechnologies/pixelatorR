#' @include generics.R
NULL

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
#' pxl_file_mpx <- minimal_mpx_pxl_file()
#' pxl_file_pna <- minimal_pna_pxl_file()
#'
#' seur_obj_mpx <- ReadMPX_Seurat(pxl_file_mpx)
#' seur_obj_mpx
#'
#' # Plot with data.frame
#' MoleculeRankPlot(seur_obj_mpx[[]])
#'
#' @export
#'
MoleculeRankPlot.data.frame <- function(
  object,
  group_by = NULL,
  ...
) {
  # Check object
  assert_non_empty_object(object, classes = "data.frame")
  if (!any(c("edges", "molecules", "n_umi") %in% colnames(object))) {
    cli::cli_abort(
      c("x" = "Either {.str edges} or {.str molecules} must be present in {.var object}")
    )
  }

  molecules_column <- intersect(colnames(object), c("molecules", "edges", "n_umi"))[1]

  assert_col_class(molecules_column, object, c("numeric", "integer"))

  if (!is.null(group_by)) {
    assert_single_value(group_by, type = "string")
    assert_col_class(group_by, object, classes = c("character", "factor"))

    if (!group_by %in% colnames(object)) {
      cli_alert_warning("Column {.str {group_by}} not found in 'object'")
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
        ggplot(., aes(rank, molecules, color = !!sym(group_by)))
      } else {
        ggplot(., aes(rank, molecules))
      }
    } +
    geom_point(size = 0.5) +
    scale_x_log10() +
    scale_y_log10() +
    labs(
      x = "Component rank (by number of molecules)",
      y = "Number of molecules"
    ) +
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
#' # Plot with Seurat object
#' MoleculeRankPlot(seur_obj_mpx)
#'
#' # Plot with Seurat object containing PNA data
#' seur_obj_pna <- ReadPNA_Seurat(pxl_file_pna)
#' MoleculeRankPlot(seur_obj_pna)
#'
#' @export
#'
MoleculeRankPlot.Seurat <- function(
  object,
  group_by = NULL,
  ...
) {
  cellrank_plot <- MoleculeRankPlot(object[[]], group_by = group_by)
  return(cellrank_plot)
}

#' Edge Rank Plot
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' The function has been replaced by \code{\link{MoleculeRankPlot}} and
#' will be removed in a future release.
#'
#' @param object A Seurat object
#' @param group_by A character specifying a column to group by
#' @param ... Additional arguments to pass to \code{\link{MoleculeRankPlot}}
#'
#' @rdname EdgeRankPlot
#' @concept plots
#'
#' @examples
#' library(pixelatorR)
#'
#' # Load example data as a Seurat object
#' pxl_file <- minimal_mpx_pxl_file()
#' seur_obj <- ReadMPX_Seurat(pxl_file)
#' seur_obj
#'
#' EdgeRankPlot(seur_obj)
#'
#' @export
#'
EdgeRankPlot <- function(
  object,
  group_by = NULL,
  ...
) {
  assert_class(object, "Seurat")

  moleculerank_plot <- MoleculeRankPlot(object[[]], group_by = group_by, ...)
  return(moleculerank_plot)
}
