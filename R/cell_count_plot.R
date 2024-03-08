#' @include generics.R
NULL

#' @rdname CellCountPlot
#' @method CellCountPlot data.frame
#'
#' @examples
#'
#' library(pixelatorR)
#' # Set arrow data output directory to temp for tests
#' options(pixelatorR.arrow_outdir = tempdir())
#'
#' # Load example data as a Seurat object
#' pxl_file <- system.file("extdata/PBMC_10_cells",
#'   "Sample01_test.pxl",
#'   package = "pixelatorR"
#' )
#' seur_obj <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)
#' seur_obj
#'
#' # Add random labels to color by
#' seur_obj$labels <- sample(c("A", "B"), size = ncol(seur_obj), replace = TRUE)
#'
#' # Plot with data.frame and color by labels
#' CellCountPlot(seur_obj[[]], color_by = "labels")
#'
#' @export
#'
CellCountPlot.data.frame <- function(
  object,
  group_by = NULL,
  color_by,
  show_count = TRUE,
  flip_axes = FALSE,
  ...
) {
  # Validate object
  if (!is.null(group_by)) {
    if (!group_by %in% colnames(object)) {
      abort(glue("'{group_by}' is missing"))
    }
    stopifnot(
      "'group_by' must be a character or factor" =
        inherits(object[, group_by, drop = TRUE],
                 what = c("character", "factor")),
      "'group_by' and 'color_by' cannot be identical" =
        group_by != color_by
    )
  }
  if (!color_by %in% colnames(object)) {
    abort(glue("'{color_by}' is missing"))
  }
  stopifnot(
    "'color_by' must be a character or factor" =
      inherits(object[, color_by, drop = TRUE],
               what = c("character", "factor"))
  )

  # Create plot
  if (!is.null(group_by)) {
    gg <- object %>%
      group_by(.data[[group_by]], .data[[color_by]]) %>%
      count() %>%
      ungroup()
    nudge_y <- max(gg$n) / 20
    p <- gg %>%
      ggplot(aes(.data[[group_by]], n, fill = .data[[color_by]])) +
      geom_col(position = position_dodge(width = 0.95)) +
      {
        if (show_count) {
          geom_text(aes(.data[[group_by]], n + nudge_y,
                        label = n,
                        group = .data[[color_by]]),
            position = position_dodge(width = 1)
          )
        }
      } +
      labs(title = paste0("Cell counts for ", group_by,
                          " colored by ", color_by))
  } else {
    gg <- object %>%
      group_by(.data[[color_by]]) %>%
      count() %>%
      ungroup()
    nudge_y <- max(gg$n) / 20
    p <- gg %>%
      ggplot(aes(.data[[color_by]], n, fill = .data[[color_by]])) +
      geom_col(position = position_dodge(width = 0.95)) +
      {
        if (show_count) {
          geom_text(aes(.data[[color_by]], n + nudge_y, label = n))
        }
      } +
      labs(title = paste0("Cell counts for ", group_by,
                          " colored by ", group_by))
  }

  # Modfy theme
  p <- p +
    scale_y_continuous(expand = expansion(c(0, 0.05))) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))

  # Handle axis flipping
  if (flip_axes) {
    p <- p + coord_flip()
  }

  return(p)
}

#' @rdname CellCountPlot
#' @method CellCountPlot Seurat
#'
#' @examples
#' # Plot with Seurat object
#' CellCountPlot(seur_obj, color_by = "labels")
#'
#' # Color by sample in merged data
#' seur_obj1 <- seur_obj2 <- seur_obj
#' seur_obj1$sample <- "1"
#' seur_obj2$sample <- "2"
#' seur_obj_merged <- merge(seur_obj1, seur_obj2, add.cell.ids = c("A", "B"))
#' CellCountPlot(seur_obj_merged, group_by = "labels", color_by = "sample")
#'
#' @export
#'
CellCountPlot.Seurat <- function (
  object,
  group_by = NULL,
  color_by,
  show_count = TRUE,
  flip_axes = FALSE,
  ...
) {
  # Extract meta.data
  mData <- object[[]]

  # Create plot
  CellCountPlot(mData,
                group_by = group_by,
                color_by = color_by,
                show_count = show_count,
                flip_axes = flip_axes)
}
