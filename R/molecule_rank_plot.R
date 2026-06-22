#' @include generics.R
NULL

#' @rdname MoleculeRankPlot
#' @method MoleculeRankPlot data.frame
#' @concept plots
#'
#' @examples
#'
#' library(pixelatorR)
#'
#' # Load example data as a Seurat object
#' pxl_file_pna <- minimal_pna_pxl_file()
#'
#' seur_obj_pna <- ReadPNA_Seurat(pxl_file_pna)
#' seur_obj_pna
#'
#' # Plot with data.frame
#' MoleculeRankPlot(seur_obj_pna[[]])
#'
#' @export
#'
MoleculeRankPlot.data.frame <- function(
  object,
  group_by = NULL,
  n_umi_min_threshold = NULL,
  n_umi_max_threshold = NULL,
  highlight_cell_counts = TRUE,
  rug = FALSE,
  split = FALSE,
  ...
) {
  # Check object
  assert_non_empty_object(object, classes = "data.frame")
  assert_single_value(n_umi_min_threshold, type = "numeric", allow_null = TRUE)
  assert_within_limits(n_umi_min_threshold, c(0, Inf), allow_null = TRUE)
  assert_single_value(n_umi_max_threshold, type = "numeric", allow_null = TRUE)
  assert_within_limits(n_umi_max_threshold, c(0, Inf), allow_null = TRUE)
  if (!is.null(n_umi_min_threshold) && !is.null(n_umi_max_threshold)) {
    if (n_umi_max_threshold <= n_umi_min_threshold) {
      cli::cli_abort(
        c("x" = "{.var n_umi_max_threshold} must be greater than {.var n_umi_min_threshold}")
      )
    }
  }
  assert_single_value(highlight_cell_counts, type = "bool")
  assert_single_value(rug, type = "bool")
  assert_single_value(split, type = "bool")
  if (!any(c("edges", "molecules", "n_umi") %in% colnames(object))) {
    cli::cli_abort(
      c("x" = "Either {.str edges} or {.str molecules} must be present in {.var object}")
    )
  }

  molecules_column <- intersect(colnames(object), c("molecules", "edges", "n_umi"))[1]

  assert_col_class(molecules_column, object, c("numeric", "integer"))

  if (!is.null(group_by)) {
    assert_single_value(group_by, type = "string")
    assert_col_in_data(group_by, object)
    assert_col_class(group_by, object, classes = c("character", "factor"))

    if (!group_by %in% colnames(object)) {
      cli_alert_warning("Column {.str {group_by}} not found in 'object'")
      group_by <- NULL
    } else {
      object <- object %>%
        mutate(!!sym(group_by) := factor(!!sym(group_by))) %>%
        group_by_at(group_by)
    }
  } else {
    if (split) {
      cli::cli_abort(
        c("x" = "Cannot split by group when {.var group_by} is NULL")
      )
    }
  }

  object <-
    object %>%
    rename(molecules = !!sym(molecules_column)) %>%
    select(molecules, any_of(group_by)) %>%
    mutate(rank = rank(-molecules, ties.method = "random"))

  # Set min/max umi thresholds if not specified
  if (is.null(n_umi_min_threshold)) {
    n_umi_min_threshold <- 0
  }
  if (is.null(n_umi_max_threshold)) {
    n_umi_max_threshold <- Inf
  }

  object <- object %>%
    mutate(
      umi_category = case_when(
        molecules < n_umi_min_threshold ~ "Low",
        molecules > n_umi_max_threshold ~ "High",
        TRUE ~ "Normal"
      )
    )

  umi_outliers <- object %>%
    {
      if (split) {
        group_by(., umi_category, !!sym(group_by))
      } else {
        group_by(., umi_category)
      }
    } %>%
    count()

  # Create edge rank plot
  cellrank_plot <- object %>%
    {
      if (!is.null(group_by)) {
        ggplot(., aes(x = rank, y = molecules, color = !!sym(group_by)))
      } else {
        ggplot(., aes(x = rank, y = molecules))
      }
    } +
    {
      if (split) {
        geom_point(size = 0.5, color = "black")
      } else {
        geom_point(size = 0.5)
      }
    } +
    scale_x_log10(limits = range(object$rank)) +
    scale_y_log10() +
    {
      # Add 1d density on left side if rug is TRUE
      if (rug) {
        geom_rug(
          linewidth = 1,
          sides = "l",
          alpha = 0.05,
          color = "black"
        )
      }
    } +
    labs(
      x = "Component rank (by number of molecules)",
      y = "Number of molecules"
    ) +
    theme_bw()

  # Add horizontal lines for min/max umi thresholds if specified
  if (n_umi_min_threshold > 0) {
    cellrank_plot <- cellrank_plot +
      geom_hline(
        yintercept = n_umi_min_threshold,
        linetype = "dashed",
        color = "red"
      )
  }
  if (is.finite(n_umi_max_threshold)) {
    cellrank_plot <- cellrank_plot +
      geom_hline(
        yintercept = n_umi_max_threshold,
        linetype = "dashed",
        color = "red"
      )
  }

  if (highlight_cell_counts) {
  
    # Define the boundary limits of the 3 zones
    min_val <- min(object$molecules, na.rm = TRUE)
    max_val <- max(object$molecules, na.rm = TRUE)

    y_min_thresh <- if (n_umi_min_threshold > 0) n_umi_min_threshold else min_val
    y_max_thresh <- if (is.finite(n_umi_max_threshold)) n_umi_max_threshold else max_val
    
    # Calculate the exact visual center for each zone (Geometric Mean)
    umi_outliers <- umi_outliers %>%
      mutate(
        y_position = case_when(
          umi_category == "Low"    ~ sqrt(as.numeric(min_val) * as.numeric(y_min_thresh)),
          umi_category == "Normal" ~ sqrt(as.numeric(y_min_thresh) * as.numeric(y_max_thresh)),
          umi_category == "High"   ~ sqrt(as.numeric(y_max_thresh) * as.numeric(max_val))
        ),
        col = case_when(
          umi_category == "Normal" ~ "#496389",
          TRUE                     ~ "#f94526"
        ),
        label_text = paste0(umi_category, ": ", n)
      )

    # Define a safe x-position (left-aligned with a bit of padding)
    fixed_x <- min(object$rank, na.rm = TRUE) * 1.5 

    cellrank_plot <- cellrank_plot +
      geom_text(
        data = umi_outliers,
        aes(
          x = fixed_x, 
          y = y_position, 
          label = label_text,
          # Dynamically handle grouping if split is TRUE
          group = if(split) !!sym(group_by) else NULL 
        ),
        color = umi_outliers$col,
        hjust = 0,             # Left-aligns text at fixed_x
        inherit.aes = FALSE    # Prevents clash with main plot aesthetics
      )
  }

  if (split) {
    cellrank_plot <- cellrank_plot +
      facet_wrap(vars(!!sym(group_by)))
  }

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
#' MoleculeRankPlot(seur_obj_pna)
#'
#' @export
#'
MoleculeRankPlot.Seurat <- function(
  object,
  group_by = NULL,
  n_umi_min_threshold = NULL,
  n_umi_max_threshold = NULL,
  highlight_cell_counts = TRUE,
  rug = FALSE,
  split = FALSE,
  ...
) {
  cellrank_plot <- MoleculeRankPlot(
    object[[]],
    group_by = group_by,
    n_umi_min_threshold = n_umi_min_threshold,
    n_umi_max_threshold = n_umi_max_threshold,
    highlight_cell_counts = highlight_cell_counts,
    rug = rug,
    split = split,
    ...
  )
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
