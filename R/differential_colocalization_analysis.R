# Declarations used in package check
globalVariables(
  names = c('pearson_z', 'p', 'p.value'),
  package = 'pixelatorR',
  add = TRUE
)
#' @include generics.R
NULL

#' @param cl  A cluster object created by makeCluster, or an integer
#' to indicate number of child-processes (integer values are ignored
#' on Windows) for parallel evaluations. See Details on performance
#' in the documentation for \code{pbapply}. The default is NULL,
#' which means that no parallelization is used.
#' @rdname RunDCA
#' @method RunDCA data.frame
#'
#' @export
#'
RunDCA.data.frame <- function (
  object,
  target,
  reference,
  contrast_column,
  group_vars = NULL,
  alternative = c("two.sided", "less", "greater"),
  conf_int = TRUE,
  p_adjust_method = c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr"),
  cl = NULL,
  verbose = TRUE,
  ...
) {

  # Validate input parameters
  stopifnot(
    "'contrast_column' must be a valid column name" =
      inherits(contrast_column, what = "character") &&
      (length(contrast_column) == 1) &&
      (contrast_column %in% colnames(object))
  )
  group_vector <- object[, contrast_column, drop = TRUE]
  stopifnot(
    "'contrast_column' must be a character vector or a factor" =
      inherits(group_vector, what = c("character", "factor"))
  )
  stopifnot(
    "'target' must be present in 'contrast_column' column" =
      inherits(target, what = "character") &&
      (length(target) == 1) &&
      (target %in% group_vector),
    "'reference' must be present in 'contrast_column' column" =
      inherits(reference, what = "character") &&
      (length(reference) == 1) &&
      (reference %in% group_vector)
  )
  stopifnot(
    "'morans_z' and 'component' must be present in colocalization score table" =
      all(c("marker_1", "marker_2", "pearson_z", "component") %in% colnames(object)),
    "'conf_int' must be TRUE or FALSE" =
      inherits(conf_int, what = "logical") &&
      (length(conf_int) == 1)
  )
  if (!is.null(group_vars)) {
    stopifnot(
      "'group_vars' must be valid column names" =
        inherits(group_vars, what = "character") &&
        (length(group_vars) >= 1) &&
        all(group_vars %in% colnames(object))
    )
    for (group_var in group_vars) {
      stopifnot(
        "'group_vars' must be character vectors or factors" =
          inherits(object[, group_var, drop = TRUE], what = c("character", "factor"))
        )
    }
  }
  alternative <- match.arg(alternative, choices = c("two.sided", "less", "greater"))
  p_adjust_method <- match.arg(p_adjust_method,
                               choices = c("bonferroni", "holm", "hochberg",
                                           "hommel", "BH", "BY", "fdr"))

  # Keep relevant columns and filter self matches
  object <- object %>%
    select(any_of(c("marker_1", "marker_2", "pearson_z", "component", contrast_column, group_vars))) %>%
    filter(marker_1 != marker_2)

  # Group data and get contrasts
  if (verbose && check_global_verbosity()) {
    cli_alert_info("Splitting data by: {paste(c('marker_1', 'marker_2', group_vars), collapse = ', ')}")
  }
  test_groups <- object %>%
    {
      if (!is.null(group_vars)) {
        group_by_at(., all_of(c("marker_1", "marker_2", group_vars)))
      } else {
        group_by_at(., c("marker_1", "marker_2"))
      }
    }

  # Get group keys
  test_groups_keys <- test_groups %>% group_keys()

  if (verbose && check_global_verbosity()) {
    cli_alert_info("Running {nrow(test_groups_keys)} tests")
  }

  # Split table by group keys
  test_groups <- test_groups %>% group_split()

  coloc_test <- pblapply(test_groups, function(colocalization_contrast) {

    # Get numeric data values
    x <- colocalization_contrast %>% filter(.data[[contrast_column]] == target) %>% pull(pearson_z)
    y <- colocalization_contrast %>% filter(.data[[contrast_column]] == reference) %>% pull(pearson_z)

    # Run wilcox.test
    result <- wilcox.test(x = x, y = y, paired = FALSE, alternative = alternative, conf.int = conf_int)

    # Tidy up results
    result <- result %>% .tidy() %>%
      mutate(n1 = length(x), n2 = length(y), method = "Wilcoxon", alternative = alternative, data_type = "pearson_z",
             target = target, reference = reference, p = signif(p.value, 3)) %>%
      select(c("estimate", "data_type", "target", "reference", "n1", "n2", "statistic",
               "p", "conf.low", "conf.high", "method", "alternative"))
    return(result)
  }, cl = cl)

  # Bind results from pol_test list
  coloc_test_bind <- do.call(bind_rows, lapply(seq_along(coloc_test), function(i) {

    # Add marker column
    coloc_test_with_groups <- coloc_test[[i]] %>% mutate(
      marker_1 = test_groups_keys[i, 1, drop = TRUE],
      marker_2 = test_groups_keys[i, 2, drop = TRUE]
    )

    # Add additional group columns
    if (!is.null(group_vars)) {
      for (group_var in group_vars) {
        coloc_test_with_groups[[group_var]] = test_groups_keys[i, group_var, drop = TRUE]
      }
    }
    return(coloc_test_with_groups)
  }))

  # Adjust p-values
  coloc_test_bind <- coloc_test_bind %>%
    mutate(p_adj = p.adjust(p, p_adjust_method))

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
#'                         "five_cells.pxl",
#'                         package = "pixelatorR")
#' # Seurat objects
#' seur1 <- seur2 <- ReadMPX_Seurat(pxl_file)
#' seur1$sample <- "Sample1"
#' seur2$sample <- "Sample2"
#' seur_merged <- merge(seur1, seur2, add.cell.ids = c("A", "B"))
#'
#' # Subset data to run test on a few markers
#' seur_merged <- subset(seur_merged,
#'                       features = c("CD3", "CD4", "CD8", "CD19",
#'                                    "CD20", "CD45RA", "CD45RO"))
#'
#' # Run DCA
#' dca_markers <- RunDCA(seur_merged, contrast_column = "sample",
#'                       target = "Sample1", reference = "Sample2")
#' dca_markers
#'
#' @export
#'
RunDCA.Seurat <- function (
  object,
  target,
  reference,
  contrast_column,
  assay = NULL,
  group_vars = NULL,
  alternative = c("two.sided", "less", "greater"),
  conf_int = TRUE,
  p_adjust_method = c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr"),
  verbose = TRUE,
  ...
) {

  # Validate input parameters
  stopifnot(
    "'contrast_column' must be a valid meta data column name" =
      inherits(contrast_column, what = "character") &&
      (length(contrast_column) == 1) &&
      (contrast_column %in% colnames(object[[]]))
  )
  group_vector <- object[[]][, contrast_column, drop = TRUE]
  stopifnot(
    "'contrast_column' must be a character vector or a factor" =
      inherits(group_vector, what = c("character", "factor"))
  )
  stopifnot(
    "'target' must be present in 'contrast_column' column" =
      inherits(target, what = "character") &&
      (length(target) == 1) &&
      (target %in% group_vector),
    "'reference' must be present in 'group_var' column" =
      inherits(reference, what = "character") &&
      (length(reference) == 1) &&
      (reference %in% group_vector)
  )
  if (!is.null(group_vars)) {
    stopifnot(
      "'group_vars' must be valid meta data column names" =
        inherits(group_vars, what = "character") &&
        (length(group_vars) >= 1) &&
        all(group_vars %in% colnames(object[[]]))
    )
    for (group_var in group_vars) {
      stopifnot(
        "'group_vars' must be character vectors or factors" =
          inherits(object[[]][, group_var, drop = TRUE], what = c("character", "factor"))
      )
    }
  }
  alternative <- match.arg(alternative, choices = c("two.sided", "less", "greater"))
  p_adjust_method <- match.arg(p_adjust_method,
                               choices = c("bonferroni", "holm", "hochberg",
                                           "hommel", "BH", "BY", "fdr"))

  # Use default assay if assay = NULL
  if (!is.null(assay)) {
    stopifnot(
      "'assay' must be a character of length 1" =
        is.character(assay) &&
        (length(assay) == 1)
    )
  } else {
    assay <- DefaultAssay(object)
  }

  # Fetch colocalization scores
  colocalization_data <- object[[assay]]@colocalization

  # Add group data
  group_data <- object[[]] %>%
    {
      if (!is.null(group_vars)) {
        select(., all_of(c(contrast_column, group_vars)))
      } else {
        select(., all_of(contrast_column))
      }
    } %>%
    as_tibble(rownames = "component")

  # Add group data to colocalization table
  colocalization_data <- colocalization_data %>%
    left_join(y = group_data, by = "component")

  # Remove redundant columns
  colocalization_data <- colocalization_data %>%
    select(-any_of(c("pearson", "pearson_mean", "pearson_stdev",
                     "pearson_p_value", "pearson_p_value_adjusted",
                     "jaccard", "jaccard_mean", "jaccard_stdev",
                     "jaccard_z", "jaccard_p_value", "jaccard_p_value_adjusted")))

  # Run DPA
  coloc_test_bind <- RunDCA(colocalization_data,
                            target = target,
                            reference = reference,
                            contrast_column = contrast_column,
                            group_vars = group_vars,
                            alternative = alternative,
                            conf_int = conf_int,
                            p_adjust_method = p_adjust_method,
                            verbose = verbose)

  return(coloc_test_bind)
}
