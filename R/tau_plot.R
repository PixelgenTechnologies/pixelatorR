# Declarations used in package check
globalVariables(
  names = c('tau', 'tau_type', 'umi_per_upia'),
  package = 'pixelatorR',
  add = TRUE
)
#' @include generics.R
NULL

#' @rdname TauPlot
#' @method TauPlot data.frame
#'
#' @examples
#'
#' library(pixelatorR)
#' # Set arrow data output directory to temp for tests
#' options(pixelatorR.arrow_outdir = tempdir())
#'
#' # Load example data as a Seurat object
#' pxl_file <- system.file("extdata/five_cells",
#'                         "five_cells.pxl",
#'                         package = "pixelatorR")
#' seur_obj <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)
#' seur_obj
#'
#' # Plot with data.frame
#' TauPlot(seur_obj[[]])
#'
#' @export
#'
TauPlot.data.frame <- function (
  object,
  group_by = NULL,
  ...
) {

  # Validate object
  if (!all(c("umi_per_upia", "tau", "tau_type") %in% colnames(object))) {
    abort("'umi_per_upia', 'tau' and 'tau_type' must be available in the data.frame")
  }
  stopifnot(
    "'umi_per_upia' must be a numeric vector" =
      is.numeric(object[, "umi_per_upia", drop = TRUE]),

    "'tau' must be a numeric vector" =
      is.numeric(object[, "tau", drop = TRUE]),

    "'tau_type' must be a factor or a character vector" =
      inherits(object[, "tau_type", drop = TRUE],
               what = c("character", "factor"))
  )
  if (!is.null(group_by)) {
    if (!group_by %in% colnames(object)) {
      abort(glue("'{group_by}' is missing"))
    }
    stopifnot("'group_by' must be a character or factor" =
                inherits(object[, group_by, drop = TRUE],
                         what = c("character", "factor")))
  }

  # Create plot
  object %>%
    ggplot(aes(tau, umi_per_upia, color = tau_type)) +
    geom_point() +
    {
      if (!is.null(group_by)) {
        facet_grid(~ sample)
      }
    } +
    scale_y_log10() +
    scale_color_manual(values = c("high" = "orangered2", "low" = "skyblue3", "normal" = "gray")) +
    theme_minimal() +
    labs(x = "Marker specificity (Tau)", y = "Pixel content (UMI/UPIA)")
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
TauPlot.Seurat <- function (
  object,
  group_by = NULL,
  ...
) {

  # Validate object
  if (!all(c("umi_per_upia", "tau", "tau_type") %in% colnames(object[[]]))) {
    abort("'umi_per_upia' and 'tau' must be available from the meta.data slot")
  }

  # Extract meta.data
  mData <- object[[]]

  # Create plot
  TauPlot(mData, group_by = group_by)
}
