#' @include generics.R
NULL

#' @rdname CellCountPlot
#' @method CellCountPlot data.frame
#'
#' @examples
#'
#' library(pixelatorR)
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
#' CellCountPlot(seur_obj[[]], color.by = "labels")
#'
#' @export
#'
CellCountPlot.data.frame <- function(
  object,
  group.by = NULL,
  color.by,
  show_count = TRUE,
  flip_axes = FALSE,
  ...
) {
  # Validate object
  if (!is.null(group.by)) {
    if (!group.by %in% colnames(object)) {
      abort(glue("'{group.by}' is missing"))
    }
    stopifnot(
      "'group.by' must be a character or factor" =
        inherits(object[, group.by, drop = TRUE],
                 what = c("character", "factor")),
      "'group.by' and 'color.by' cannot be identical" =
        group.by != color.by
    )
  }
  if (!color.by %in% colnames(object)) {
    abort(glue("'{color.by}' is missing"))
  }
  stopifnot(
    "'color.by' must be a character or factor" =
      inherits(object[, color.by, drop = TRUE],
               what = c("character", "factor"))
  )

  # Create plot
  if (!is.null(group.by)) {
    gg <- object %>%
      group_by(.data[[group.by]], .data[[color.by]]) %>%
      count() %>%
      ungroup()
    nudge_y <- max(gg$n) / 20
    p <- gg %>%
      ggplot(aes(.data[[group.by]], n, fill = .data[[color.by]])) +
      geom_col(position = position_dodge(width = 0.95)) +
      {
        if (show_count) {
          geom_text(aes(.data[[group.by]], n + nudge_y,
                        label = n,
                        group = .data[[color.by]]),
            position = position_dodge(width = 1)
          )
        }
      } +
      labs(title = paste0("Cell counts for ", group.by,
                          " colored by ", color.by))
  } else {
    gg <- object %>%
      group_by(.data[[color.by]]) %>%
      count() %>%
      ungroup()
    nudge_y <- max(gg$n) / 20
    p <- gg %>%
      ggplot(aes(.data[[color.by]], n, fill = .data[[color.by]])) +
      geom_col(position = position_dodge(width = 0.95)) +
      {
        if (show_count) {
          geom_text(aes(.data[[color.by]], n + nudge_y, label = n))
        }
      } +
      labs(title = paste0("Cell counts for ", group.by,
                          " colored by ", group.by))
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
#' CellCountPlot(seur_obj, color.by = "labels")
#'
#' # Color by sample in merged data
#' seur_obj1 <- seur_obj2 <- seur_obj
#' seur_obj1$sample <- "1"
#' seur_obj2$sample <- "2"
#' seur_obj_merged <- merge(seur_obj1, seur_obj2, add.cell.ids = c("A", "B"))
#' CellCountPlot(seur_obj_merged, group.by = "labels", color.by = "sample")
#'
#' @export
#'
CellCountPlot.Seurat <- function (
  object,
  group.by = NULL,
  color.by,
  show_count = TRUE,
  flip_axes = FALSE,
  ...
) {
  # Extract meta.data
  mData <- object[[]]

  # Create plot
  CellCountPlot(mData,
                group.by = group.by,
                color.by = color.by,
                show_count = show_count,
                flip_axes = flip_axes)
}
