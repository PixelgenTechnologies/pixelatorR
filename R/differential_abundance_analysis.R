#' @include generics.R
#' @include utils.R
#' @include differential_analysis_utils.R
NULL

#' @param assay Name of assay to use
#' @param mean_fxn Function to use for fold change or average difference calculation
#' in \code{FindMarkers}. See documentation for the \code{mean.fxn} parameter in
#' \code{FindMarkers} for more information.
#' @param fc_name Name of the fold change, average difference, or custom function
#' column in the output tibble. See documentation for the \code{fc.name} parameter in
#' \code{FindMarkers} for more information.
#'
#' @rdname RunDAA
#' @method RunDAA Seurat
#'
#' @examples
#' library(pixelatorR)
#' library(dplyr)
#' library(SeuratObject)
#'
#' pxl_file <- system.file("extdata/five_cells",
#'   "five_cells.pxl",
#'   package = "pixelatorR"
#' )
#' # Seurat objects
#' se <- ReadMPX_Seurat(pxl_file)
#' se <- merge(se, rep(list(se), 9), add.cell.ids = LETTERS[1:10])
#' se$sample <- c("T", "C", "C", "C", "C") %>% rep(times = 10)
#' se <- Seurat::NormalizeData(se, normalization.method = "CLR", margin = 2)
#'
#' # Run DPA
#' daa_markers <- RunDAA(se,
#'   contrast_column = "sample",
#'   targets = "T", reference = "C"
#' )
#' daa_markers
#'
#' @export
#'
RunDAA.Seurat <- function(
  object,
  contrast_column,
  reference,
  targets = NULL,
  assay = NULL,
  group_vars = NULL,
  mean_fxn = rowMeans,
  fc_name = "difference",
  p_adjust_method = c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr"),
  verbose = TRUE,
  ...
) {
  expect_Seurat()

  # Validate input parameters
  stopifnot(
    "'contrast_column' must be available in Seurat object meta.data" =
      contrast_column %in% colnames(object[[]])
  )

  .validate_da_input(object[[]], contrast_column, reference, targets, group_vars)

  # Use default assay if assay = NULL
  assay <- .validate_or_set_assay(object, assay)

  # Define targets as all groups except the reference if not specified
  targets <- targets %||% setdiff(unique(object[[]][, contrast_column, drop = TRUE]), reference)

  # Check multiple choice args
  p_adjust_method <- match.arg(p_adjust_method,
    choices = c(
      "bonferroni", "holm", "hochberg",
      "hommel", "BH", "BY", "fdr"
    )
  )

  # Add group data
  group_data <- .select_group_data(object[[]], contrast_column, group_vars)

  # Group data and get contrasts
  if (verbose && check_global_verbosity() && !is.null(group_vars)) {
    cli_alert_info("Splitting data by: {paste(c(group_vars), collapse = ', ')}")
  }

  test_groups <- group_data %>%
    {
      if (!is.null(group_vars)) {
        group_by(., pick(all_of(group_vars)))
      } else {
        .
      }
    }

  # Get group keys
  test_groups_keys <- test_groups %>% group_keys()
  if (nrow(test_groups_keys) == 1) {
    group_data_split <- list(test_groups)
  } else {
    group_data_split <- test_groups %>%
      group_split()
  }

  if (verbose && check_global_verbosity() && !is.null(group_vars)) {
    .print_da_group_strategy(test_groups_keys, targets, reference)
  }

  # Get assay and convert if necessary
  cg_assay <- object[[assay]]
  if (inherits(cg_assay, "CellGraphAssay")) {
    cg_assay <- as(cg_assay, "Assay")
  }
  if (inherits(cg_assay, "CellGraphAssay5")) {
    cg_assay <- as(cg_assay, "Assay5")
  }

  de_results_all <- lapply(seq_along(group_data_split), function(i) {
    reference_cells <- group_data_split[[i]] %>%
      filter(!!sym(contrast_column) == reference) %>%
      pull(component)

    # Iterate over targets
    de_results_targets <- lapply(targets, function(target) {
      target_cells <- group_data_split[[i]] %>%
        filter(!!sym(contrast_column) == target) %>%
        pull(component)

      cur_assay_subset <- subset(cg_assay, cells = group_data_split[[i]]$component)
      target_cells <- group_data_split[[i]] %>%
        filter(!!sym(contrast_column) == target) %>%
        pull(component)
      de_results_cur <- Seurat::FindMarkers(
        cur_assay_subset,
        cells.1 = target_cells,
        cells.2 = reference_cells,
        mean.fxn = mean_fxn,
        fc.name = fc_name,
        logfc.threshold = 0,
        ...
      )

      # Add marker columns
      de_results_cur <- de_results_cur %>%
        mutate(marker = rownames(.)) %>%
        relocate(marker, .before = "p_val") %>%
        as_tibble() %>%
        mutate(target = target, reference = reference)

      # Add additional group columns
      if (!is.null(group_vars)) {
        for (group_var in group_vars) {
          de_results_cur[[group_var]] <- test_groups_keys[i, group_var, drop = TRUE]
        }
      }

      return(de_results_cur)
    }) %>%
      bind_rows()

    return(de_results_targets)
  }) %>%
    bind_rows()

  # Adjust p-values
  de_results_all <- de_results_all %>%
    rename(p = p_val, pct_1 = pct.1, pct_2 = pct.2) %>%
    select(-p_val_adj) %>%
    mutate(p_adj = p.adjust(p, p_adjust_method)) %>%
    relocate(p_adj, .after = "p")

  return(de_results_all)
}
