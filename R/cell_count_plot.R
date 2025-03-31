#' @include generics.R
NULL

#' @rdname CellCountPlot
#' @method CellCountPlot data.frame
#' @examples
#'
#' library(pixelatorR)
#'
#' # Load example data as a Seurat object
#' pxl_file <- minimal_mpx_pxl_file()
#' seur_obj <- ReadMPX_Seurat(pxl_file)
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
  as_frequency = FALSE,
  stack = FALSE,
  ...
) {
  # Validate input
  assert_col_in_data(color_by, object)
  assert_col_class(color_by, object, classes = c("character", "factor"))
  assert_single_value(group_by, type = "string", allow_null = TRUE)
  assert_single_values_are_different(group_by, color_by, allow_null = TRUE)
  assert_col_in_data(group_by, object, allow_null = TRUE)
  assert_col_class(group_by, object, classes = c("character", "factor"), allow_null = TRUE)

  # Create plot
  if (!is.null(group_by)) {
    gg <- object %>%
      group_by(.data[[group_by]], .data[[color_by]]) %>%
      count() %>%
      ungroup() %>%
      group_by(!!sym(group_by)) %>%
      mutate(frequency = n / sum(n) * 100) %>%
      ungroup()

    nudge_y <- if (as_frequency) max(gg$frequency) / 20 else max(gg$n) / 20

    p <- ggplot(gg, aes(
      x = .data[[group_by]],
      y = if (as_frequency) frequency else n,
      fill = .data[[color_by]]
    )) +
      geom_col(position = if (stack) "stack" else position_dodge(width = 0.95)) +
      {
        if (show_count) {
          geom_text(
            aes(
              label = if (as_frequency) sprintf("%.1f%%", frequency) else n,
              group = .data[[color_by]]
            ),
            position = if (stack) position_stack(vjust = 0.5) else position_dodge(width = 0.95),
            size = 3
          )
        }
      } +
      labs(
        title = ifelse(as_frequency,
          paste0("Cell frequencies for ", group_by, " colored by ", color_by),
          paste0("Cell counts for ", group_by, " colored by ", color_by)
        ),
        y = ifelse(as_frequency, "Frequency (%)", "Count")
      )
  } else {
    gg <- object %>%
      group_by(.data[[color_by]]) %>%
      count() %>%
      ungroup() %>%
      mutate(frequency = n / sum(n) * 100)

    nudge_y <- if (as_frequency) max(gg$frequency) / 20 else max(gg$n) / 20

    p <- ggplot(gg, aes(
      x = .data[[color_by]],
      y = if (as_frequency) frequency else n,
      fill = .data[[color_by]]
    )) +
      geom_col(position = position_dodge(width = 0.95)) +
      {
        if (show_count) {
          geom_text(aes(label = if (as_frequency) sprintf("%.1f%%", frequency) else n),
            position = position_dodge(width = 0.95),
            vjust = -0.5,
            size = 3
          )
        }
      } +
      labs(
        title = ifelse(as_frequency,
          paste0("Cell frequencies for ", color_by),
          paste0("Cell counts for ", color_by)
        ),
        y = ifelse(as_frequency, "Frequency (%)", "Count")
      )
  }

  # Modify theme
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
CellCountPlot.Seurat <- function(
  object,
  group_by = NULL,
  color_by,
  show_count = TRUE,
  flip_axes = FALSE,
  as_frequency = FALSE,
  stack = FALSE,
  ...
) {
  # Extract meta.data
  mData <- object[[]]

  # Create plot
  CellCountPlot(mData,
    group_by = group_by,
    color_by = color_by,
    show_count = show_count,
    flip_axes = flip_axes,
    as_frequency = as_frequency,
    stack = stack
  )
}
