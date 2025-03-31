#' @include generics.R
NULL

#' @rdname TauPlot
#' @method TauPlot data.frame
#'
#' @examples
#'
#' library(pixelatorR)
#'
#' # Load example data as a Seurat object
#' pxl_file <- minimal_mpx_pxl_file()
#' seur_obj <- ReadMPX_Seurat(pxl_file)
#' seur_obj
#'
#' # Plot with data.frame
#' TauPlot(seur_obj[[]])
#'
#' @export
#'
TauPlot.data.frame <- function(
  object,
  group_by = NULL,
  ...
) {
  # Validate object
  mol_per_upia <- intersect(c("umi_per_upia", "mean_molecules_per_a_pixel"), colnames(object))
  if (length(mol_per_upia) == 0) {
    cli::cli_abort(
      c(
        "i" = "Either {.str umi_per_upia} or {.str mean_molecules_per_a_pixel}",
        " " = "must be available in the {.cls {class(object)}}",
        "x" = "Missing columns in the {.cls {class(object)}}"
      )
    )
  }
  if (length(mol_per_upia) > 1) {
    mol_per_upia <- mol_per_upia[1]
  }
  if (!is.numeric(object[, mol_per_upia, drop = TRUE])) {
    abort(glue("'{mol_per_upia}' must be a numeric vector"))
  }
  assert_col_in_data("tau", object)
  assert_col_in_data("tau_type", object)
  assert_col_class("tau", object, classes = "numeric")
  assert_col_class("tau_type", object, classes = c("character", "factor"))
  if (!is.null(group_by)) {
    assert_col_in_data(group_by, object)
    assert_col_class(group_by, object, classes = c("character", "factor"))
  }

  y_lab <- switch(mol_per_upia,
    umi_per_upia = "Pixel content (UMI / UPIA)",
    mean_molecules_per_a_pixel = "Pixel content (mean molecules / UPIA)"
  )

  # Create plot
  object %>%
    ggplot(aes(tau, !!sym(mol_per_upia), color = tau_type)) +
    geom_point() +
    {
      if (!is.null(group_by)) {
        facet_wrap(as.formula(paste("~", group_by)))
      }
    } +
    scale_y_log10() +
    scale_color_manual(values = c("high" = "orangered2", "low" = "skyblue3", "normal" = "gray")) +
    theme_minimal() +
    labs(x = "Marker specificity (Tau)", y = y_lab)
}

#' @rdname TauPlot
#' @method TauPlot Seurat
#'
#' @examples
#' # Plot with Seurat object
#' TauPlot(seur_obj)
#'
#' # Group by sample in merged data
#' seur_obj1 <- seur_obj2 <- seur_obj
#' seur_obj1$sample <- "1"
#' seur_obj2$sample <- "2"
#' seur_obj_merged <- merge(seur_obj1, seur_obj2, add.cell.ids = c("A", "B"))
#' TauPlot(seur_obj_merged, group_by = "sample")
#'
#' @export
#'
TauPlot.Seurat <- function(
  object,
  group_by = NULL,
  ...
) {
  # Validate object
  mol_per_upia <- intersect(c("umi_per_upia", "mean_molecules_per_a_pixel"), colnames(object[[]]))
  if (length(mol_per_upia) == 0) {
    cli::cli_abort(
      c(
        "i" = "Either {.str umi_per_upia} or {.str mean_molecules_per_a_pixel}",
        " " = "must be available in the {.cls {class(object)}}",
        "x" = "Missing columns in the {.cls {class(object)}}"
      )
    )
  }
  assert_col_in_data("tau", object[[]])
  assert_col_in_data("tau_type", object[[]])

  # Extract meta.data
  mData <- object[[]]

  # Create plot
  TauPlot(mData, group_by = group_by)
}
