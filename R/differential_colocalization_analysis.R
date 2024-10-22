#' @include generics.R
#' @include utils.R
#' @include differential_analysis_utils.R
NULL

#' @param cl A cluster object created by \code{\link[parallel]{makeCluster}},
#' or an integer to indicate number of child-processes (integer values are
#' ignored on Windows) for parallel evaluations. See Details on performance
#' in the documentation for \code{pbapply}. The default is NULL,
#' which means that no parallelization is used.
#' Note that warnings are not caught when using parallel processing.
#' @rdname RunDCA
#' @method RunDCA data.frame
#'
#' @export
#'
RunDCA.data.frame <- function(
  object,
  contrast_column,
  reference,
  targets = NULL,
  group_vars = NULL,
  coloc_metric = "pearson_z",
  min_n_obs = 0,
  alternative = c("two.sided", "less", "greater"),
  conf_int = TRUE,
  p_adjust_method = c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr"),
  cl = NULL,
  verbose = TRUE,
  ...
) {
  # Validate input parameters
  .validate_da_input(
    object, contrast_column, reference, targets, group_vars,
    coloc_metric, min_n_obs, conf_int, cl,
    data_type = "colocalization"
  )

  # Define targets as all groups except the reference if not specified
  targets <- targets %||% setdiff(unique(object[, contrast_column, drop = TRUE]), reference)

  # Check multiple choice args
  alternative <- match.arg(alternative, choices = c("two.sided", "less", "greater"))
  p_adjust_method <- match.arg(p_adjust_method,
    choices = c(
      "bonferroni", "holm", "hochberg",
      "hommel", "BH", "BY", "fdr"
    )
  )

  # Keep relevant columns and filter self matches
  object <- object %>%
    select(any_of(c("marker_1", "marker_2", coloc_metric, "component", contrast_column, group_vars))) %>%
    filter(marker_1 != marker_2)

  # Group data and get contrasts
  if (verbose && check_global_verbosity()) {
    cli_alert_info("Splitting data by: {paste(c('marker_1', 'marker_2', group_vars), collapse = ', ')}")
  }
  test_groups <- object %>%
    {
      if (!is.null(group_vars)) {
        group_by(., pick(all_of(c("marker_1", "marker_2", group_vars))))
      } else {
        group_by(., pick(all_of(c("marker_1", "marker_2"))))
      }
    }

  if (min_n_obs > 0) {
    # Filter test groups by minimum number of observations allowed
    test_groups <- test_groups %>%
      group_by(!!sym(contrast_column), .add = TRUE) %>%
      mutate(n = n()) %>%
      filter(n > min_n_obs) %>%
      ungroup(!!sym(contrast_column)) %>%
      mutate(
        ref_n = sum(!!sym(contrast_column) == reference),
        target_n = sum(!!sym(contrast_column) %in% targets)
      ) %>%
      filter(ref_n > 0, target_n > 0)
    if (nrow(test_groups) == 0) {
      abort(glue(
        "Found no groups with at least {min_n_obs} observations."
      ))
    }
  }

  # Get group keys
  test_groups_keys <- test_groups %>% group_keys()

  if (verbose && check_global_verbosity()) {
    .print_da_group_strategy(test_groups_keys, targets, reference)
  }

  # Split table by group keys
  test_groups <- test_groups %>% group_split()

  # Chunk the data to enable parallel processing
  # Decide how to split the data depending on the cl
  if (inherits(cl, "cluster")) {
    # Load dplyr in each cluster
    clusterEvalQ(cl, {
      library(dplyr)
      library(rlang)
      library(glue)
      library(cli)
    })
    # Export variables to each cluster
    clusterExport(cl, c(
      "contrast_column", "targets", "reference",
      "coloc_metric", "evaluate_with_catch",
      "alternative", "conf_int", ".tidy", "group_vars"
    ),
    envir = current_env()
    )
    # Parallel processing on cluster
    if ((length(test_groups) * length(cl)) > (length(cl) * 100)) {
      # Cut into even chunks of 100 values in each chunk
      chunks <- ceiling(seq_along(test_groups) / 100)
    } else {
      chunks <- cut(seq_along(test_groups), length(cl))
    }
  } else if (is.numeric(cl)) {
    # Use sequential processing when cl is 1
    if (cl == 1) {
      chunks <- seq_along(test_groups)
    } else {
      # Parallel processing when cl is the number of child processes
      if ((length(test_groups) * cl) > (cl * 100)) {
        # Cut into even chunks of 100 values in each chunk
        chunks <- ceiling(seq_along(test_groups) / 100)
      } else {
        chunks <- cut(seq_along(test_groups), cl)
      }
    }
  } else {
    # Sequential processing when cl is NULL
    chunks <- seq_along(test_groups)
  }

  # Split test_groups into chunks
  test_groups_chunked <- split(test_groups, chunks)

  # Process chunks
  coloc_test <- pblapply(test_groups_chunked, function(test_groups_chunk) {
    coloc_test_chunk <- lapply(test_groups_chunk, function(colocalization_contrast) {
      # Fetch marker for current comparison
      marker_1 <- colocalization_contrast$marker_1 %>% unique()
      marker_2 <- colocalization_contrast$marker_2 %>% unique()

      # Get numeric data values
      x_list <- lapply(targets, function(target) {
        colocalization_contrast %>%
          filter(.data[[contrast_column]] == target) %>%
          pull(all_of(coloc_metric))
      }) %>%
        set_names(nm = targets)
      y <- colocalization_contrast %>%
        filter(.data[[contrast_column]] == reference) %>%
        pull(all_of(coloc_metric))

      # Run wilcox.test for all targets vs reference
      results <- lapply(names(x_list), function(target) {
        x <- x_list[[target]]
        # Run wilcox.test
        result <- evaluate_with_catch({
          wilcox_res <- wilcox.test(
            x = x,
            y = y,
            paired = FALSE,
            alternative = alternative,
            conf.int = conf_int
          )
        })

        if (!is.null(result$error)) {
          warn(glue("Failed to compute Wilcoxon test for marker '{marker_1}/{marker_2}': {target} vs {reference}\n",
            "  Got the following error message when running wilcox.test:\n",
            "  {col_red(result$error)}\n",
            "  This test will be skipped.",
            .trim = FALSE
          ))
          return(NULL)
        }
        if (!is.null(result$warning)) {
          warn(glue("Got the following message when running ",
            "wilcox.test test for marker '{marker_1}/{marker_2}': {target} vs {reference}\n",
            "  {col_red(result$warning)}\n",
            .trim = FALSE
          ))
        }
        result <- result$result

        # Tidy up results
        result <- result %>%
          .tidy() %>%
          mutate(
            n1 = length(x), n2 = length(y), method = "Wilcoxon",
            alternative = alternative, data_type = "pearson_z",
            target = target, reference = reference, p = signif(p.value, 3)
          ) %>%
          select(c(
            "estimate", "data_type", "target", "reference",
            "n1", "n2", "statistic", "p", "conf.low", "conf.high",
            "method", "alternative"
          )) %>%
          mutate(marker_1 = marker_1, marker_2 = marker_2)

        # Add additional group columns
        if (!is.null(group_vars)) {
          for (group_var in group_vars) {
            result[[group_var]] <- colocalization_contrast[, group_var, drop = TRUE] %>% unique()
          }
        }

        return(result)
      }) %>%
        do.call(bind_rows, .)
    })

    return(coloc_test_chunk %>% do.call(bind_rows, .))
  }, cl = cl)

  # Bind results from coloc_test list
  coloc_test_bind <- do.call(bind_rows, coloc_test)

  # Adjust p-values
  coloc_test_bind <- coloc_test_bind %>%
    mutate(p_adj = p.adjust(p, p_adjust_method)) %>%
    relocate(p_adj, .after = "p")

  return(coloc_test_bind)
}


#' @param assay Name of assay to use
#'
#' @rdname RunDCA
#' @method RunDCA Seurat
#'
#' @examples
#' library(pixelatorR)
#' library(dplyr)
#'
#' pxl_file <- system.file("extdata/five_cells",
#'   "five_cells.pxl",
#'   package = "pixelatorR"
#' )
#' # Seurat objects
#' seur1 <- seur2 <- ReadMPX_Seurat(pxl_file)
#' seur1$sample <- "Sample1"
#' seur2$sample <- "Sample2"
#' seur_merged <- merge(seur1, seur2, add.cell.ids = c("A", "B"))
#'
#' # Subset data to run test on a few markers
#' seur_merged <- subset(seur_merged,
#'   features = c(
#'     "CD3", "CD4", "CD8", "CD19",
#'     "CD20", "CD45RA", "CD45RO"
#'   )
#' )
#'
#' # Run DCA
#' dca_markers <- RunDCA(seur_merged,
#'   contrast_column = "sample",
#'   target = "Sample1", reference = "Sample2"
#' )
#' dca_markers
#'
#' @export
#'
RunDCA.Seurat <- function(
  object,
  contrast_column,
  reference,
  targets = NULL,
  assay = NULL,
  group_vars = NULL,
  coloc_metric = "pearson_z",
  min_n_obs = 0,
  alternative = c("two.sided", "less", "greater"),
  conf_int = TRUE,
  p_adjust_method = c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr"),
  cl = NULL,
  verbose = TRUE,
  ...
) {
  # Validate input parameters
  stopifnot(
    "'contrast_column' must be available in Seurat object meta.data" =
      contrast_column %in% colnames(object[[]])
  )

  # Use default assay if assay = NULL
  assay <- .validate_or_set_assay(object, assay)

  # Fetch colocalization scores
  colocalization_data <- object[[assay]]@colocalization

  # Add group data
  group_data <- .select_group_data(object[[]], contrast_column, group_vars)

  # Add group data to colocalization table
  colocalization_data <- colocalization_data %>%
    left_join(y = group_data, by = "component")

  # Run DCA
  coloc_test_bind <- RunDCA(colocalization_data,
    targets = targets,
    reference = reference,
    contrast_column = contrast_column,
    group_vars = group_vars,
    coloc_metric = coloc_metric,
    min_n_obs = min_n_obs,
    alternative = alternative,
    conf_int = conf_int,
    p_adjust_method = p_adjust_method,
    cl = cl,
    verbose = verbose
  )

  return(coloc_test_bind)
}
