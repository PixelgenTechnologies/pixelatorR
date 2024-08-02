# Declarations used in package check
globalVariables(
  names = c('morans_z', 'p', 'p.value'),
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
#' @rdname RunDPA
#' @method RunDPA data.frame
#'
#' @examples
#' library(pixelatorR)
#' library(dplyr)
#'
#' pxl_file <- system.file("extdata/five_cells",
#'                         "five_cells.pxl",
#'                         package = "pixelatorR")
#'
#' # Load polarization scores
#' polarization_table1 <- polarization_table2 <- ReadMPX_polarization(pxl_file)
#' polarization_table1$sample <- "Sample1"
#' polarization_table2$sample <- "Sample2"
#' polarization_table_merged <-  bind_rows(polarization_table1, polarization_table2)
#'
#' # Run DPA using table as input
#' dpa_markers <- RunDPA(polarization_table_merged, contrast_column = "sample",
#'                       target = "Sample1", reference = "Sample2")
#' dpa_markers
#'
#' @export
#'
RunDPA.data.frame <- function (
  object,
  contrast_column,
  reference,
  targets = NULL,
  group_vars = NULL,
  polarity_metric = c("morans_z", "morans_i"),
  min_n_obs = 0,
  cl = NULL,
  alternative = c("two.sided", "less", "greater"),
  conf_int = TRUE,
  p_adjust_method = c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr"),
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
  polarity_metric <- match.arg(polarity_metric, choices = c("morans_z", "morans_i"))
  stopifnot(
    "'polarity_metric' must be present in polarization score table" =
      polarity_metric %in% colnames(object)
  )

  group_vector <- object[, contrast_column, drop = TRUE]

  stopifnot(
    "'contrast_column' must be a character vector or a factor" =
      inherits(group_vector, what = c("character", "factor")),
    "'contrast_column' must have at least 2 groups" =
      length(unique(group_vector)) > 1
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

  if (!is.null(targets)) {
    stopifnot(
      "'targets' must be present in 'contrast_column' column" =
        inherits(targets, what = "character") &&
        all(targets %in% group_vector)
    )
  }
  targets <- targets %||% setdiff(unique(group_vector), reference)

  stopifnot(
    "'morans_z' and 'component' must be present in polarization score table" =
      all(c("marker", "morans_z", "component") %in% colnames(object)),
    "'polarity_metric' and 'component' must be present in polarization score table" =
      all(c("marker", polarity_metric, "component") %in% colnames(object)),
    "'conf_int' must be TRUE or FALSE" =
      inherits(conf_int, what = "logical") & (length(conf_int) == 1)
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
          inherits(object[, group_var, drop = TRUE],
                   what = c("character", "factor"))
      )
    }
  }

  # Check multiple choice args
  alternative <- match.arg(alternative, choices = c("two.sided", "less", "greater"))
  p_adjust_method <- match.arg(p_adjust_method,
                               choices = c("bonferroni", "holm", "hochberg",
                                           "hommel", "BH", "BY", "fdr"))

  # Discard unused data columns
  object <- object %>%
    select(all_of(c("marker", polarity_metric, "component", contrast_column, group_vars)))

  # Group data and get contrasts
  if (verbose && check_global_verbosity()) {
    cli_alert_info("Splitting data by: {paste(c('marker', group_vars), collapse = ', ')}")
    cli_alert_info("Polarity metric: '{polarity_metric}'")
  }
  test_groups <- object %>%
    {
      if (!is.null(group_vars)) {
        group_by_at(., all_of(c("marker", group_vars)))
      } else {
        group_by(., .data[["marker"]])
      }
    }

  if (min_n_obs > 0) {
    # Filter test groups by minimum number of observations allowed
    test_groups <- test_groups %>%
      group_by(!! sym(contrast_column), .add = TRUE) %>%
      mutate(n = n()) %>%
      filter(n > min_n_obs) %>%
      ungroup(!! sym(contrast_column)) %>%
      mutate(ref_n = sum(!! sym(contrast_column) == reference),
             target_n = sum(!! sym(contrast_column) %in% targets)) %>%
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
    if (length(targets) > 1) {
      cli_alert_info(glue(
        "Running {nrow(test_groups_keys)} tests for each ",
        "of the following comparisons:\n"
      ))
    } else {
      cli_alert_info(glue(
        "Running {nrow(test_groups_keys)} tests for the following comparison:\n"
      ))
    }
    print(glue("  - {paste(targets, reference, sep = ' vs ')}"))
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
      "test_groups", "test_groups_keys", "contrast_column",
      "target", "reference", "alternative", "conf_int",
      ".tidy", "group_vars", "targets", "evaluate_with_catch"
    ),
    envir = current_env()
    )
    chunks <- split(seq_along(test_groups), cut(seq_along(test_groups), length(cl)))
  } else if (is.numeric(cl)) {
    # Parallel processing
    chunks <- split(seq_along(test_groups), cut(seq_along(test_groups), cl))
  } else {
    # Sequential processing
    chunks <- as.list(seq_along(test_groups))
    cl <- NULL
  }

  # Process chunks
  pol_test <- pbapply::pblapply(chunks, function(chunk) {

    test_groups_chunk <- test_groups[chunk]
    test_groups_keys_chunk <- test_groups_keys[chunk, ]

    # Process current chunk
    pol_test_chunk <- lapply(seq_along(test_groups_chunk), function(i) {

      polarity_contrast <- test_groups_chunk[[i]]

      # Fetch marker for current comparison
      marker <- test_groups_keys_chunk[i, "marker", drop = TRUE]

      # Get numeric data values for all targets
      x_list <- lapply(targets, function(target) {
        polarity_contrast %>%
          filter(.data[[contrast_column]] == target) %>%
          pull(all_of(polarity_metric))
      }) %>% set_names(nm = targets)
      y <- polarity_contrast %>%
        filter(.data[[contrast_column]] == reference) %>%
        pull(all_of(polarity_metric))

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
          warn(glue("Failed to compute Wilcoxon test for marker '{marker}': {target} vs {reference}\n",
                    "  Got the following error message when running wilcox.test:\n",
                    "  {col_red(result$error)}\n",
                    "  This test will be skipped.", .trim = FALSE))
          return(NULL)
        }
        if (!is.null(result$warning)) {
          warn(glue("Got the following message when running ",
                    "wilcox.test test for marker '{marker}': {target} vs {reference}\n",
                    "  {col_red(result$warning)}\n", .trim = FALSE))
        }
        result <- result$result

        # Tidy up results
        result <- result %>%
          .tidy() %>%
          mutate(
            n1 = length(x), n2 = length(y), method = "Wilcoxon",
            alternative = alternative, data_type = polarity_metric,
            target = target, reference = reference, p = signif(p.value, 3)
          ) %>%
          select(c(
            "estimate", "data_type", "target", "reference", "n1", "n2", "statistic",
            "p", "conf.low", "conf.high", "method", "alternative"
          )) %>%
          mutate(marker = marker)

        # Add additional group columns
        if (!is.null(group_vars)) {
          for (group_var in group_vars) {
            result[[group_var]] <- test_groups_keys[i, group_var, drop = TRUE]
          }
        }

        return(result)
      }) %>% do.call(bind_rows, .)

    })

    return(pol_test_chunk %>% do.call(bind_rows, .))
  }, cl = cl)

  pol_test_bind <- do.call(bind_rows, pol_test)

  # Adjust p-values
  pol_test_bind <- pol_test_bind %>%
    mutate(p_adj = p.adjust(p, p_adjust_method))

  return(pol_test_bind)
}


#' @param assay Name of assay to use
#'
#' @rdname RunDPA
#' @method RunDPA Seurat
#'
#' @examples
#' # Seurat objects
#' seur1 <- seur2 <- ReadMPX_Seurat(pxl_file)
#' seur1$sample <- "Sample1"
#' seur2$sample <- "Sample2"
#' seur_merged <- merge(seur1, seur2, add.cell.ids = c("A", "B"))
#'
#' # Run DPA
#' dpa_markers <- RunDPA(seur_merged, contrast_column = "sample",
#'                       target = "Sample1", reference = "Sample2")
#' dpa_markers
#'
#' @export
#'
RunDPA.Seurat <- function (
  object,
  contrast_column,
  reference,
  targets = NULL,
  assay = NULL,
  group_vars = NULL,
  polarity_metric = c("morans_z", "morans_i"),
  min_n_obs = 0,
  cl = NULL,
  alternative = c("two.sided", "less", "greater"),
  conf_int = TRUE,
  p_adjust_method = c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr"),
  verbose = TRUE,
  ...
) {

  # Validate input parameters
  stopifnot(
    "'contrast_column' must be available in Seurat object meta.data" =
      contrast_column %in% colnames(object[[]])
  )
  polarity_metric <- match.arg(polarity_metric, choices = c("morans_z", "morans_i"))

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

  # Fetch polarization scores
  polarization_data <- object[[assay]]@polarization

  # Add group data
  group_data <- object[[]] %>%
    {
      if (!is.null(group_vars)) {
        stopifnot(
          "'group_vars' must be valid meta data column names" =
            inherits(group_vars, what = "character") &&
              (length(group_vars) >= 1) &&
              all(group_vars %in% colnames(.))
        )
        select(., all_of(c(contrast_column, group_vars)))
      } else {
        select(., all_of(contrast_column))
      }
    } %>%
    as_tibble(rownames = "component")

  # Add group data to polarization table
  polarization_data <- polarization_data %>%
    left_join(y = group_data, by = "component")

  # Remove redundant columns
  polarization_data <- polarization_data %>%
    select(-any_of(c("morans_p_value", "morans_p_adjusted",
                     setdiff(c("morans_z", "morans_i"), polarity_metric))))

  # Run DPA
  pol_test_bind <- RunDPA(polarization_data,
    targets = targets,
    reference = reference,
    contrast_column = contrast_column,
    group_vars = group_vars,
    polarity_metric = polarity_metric,
    min_n_obs = min_n_obs,
    cl = cl,
    alternative = alternative,
    conf_int = conf_int,
    p_adjust_method = p_adjust_method,
    verbose = verbose
  )

  return(pol_test_bind)
}
