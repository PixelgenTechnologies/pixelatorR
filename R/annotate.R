#' Automatic annotation of cell types
#'
#' This function finds anchors in a reference dataset and transfers the cell type
#' of the reference to a query Seurat dataset. Two columns of cell type annotations
#' are added to the metadata of the query dataset, one high-level and one more
#' fine-grained. In addition, one can optionally summarize these annotations
#' for a specific clustering solution. When choosing this, the most common annotation
#' per cluster is added to the metadata of the object.
#'
#' @param object Seurat object to which you want to add celltype annotation
#' @param reference Seurat object that contains annotated reference data, see examples
#' for a programmatic way to obtain a PBMC reference dataset
#' @param summarize_by_column If NULL, return only the cell type annotation per cell.
#' If the name of a column in the object metadata is provided, the annotations will
#' be summarized to the factors of this column, i.e. the most common annotation per
#' factor level will be retained. Useful for annotating the output of the Seurat
#' clustering, for example.
#' @param reference_assay Assay in the reference dataset used to find the anchors,
#' Default: 'ADT'
#' @param query_assay Assay in the query dataset used to find the anchors. It should
#' probably be set to a normalized assay that was run on the object (e.g. "dsb"),
#' Default: 'PNA'
#' @param reference_groups A character vector of reference groups to use for annotation.
#' @param normalization_method normalization method used during \code{FindTransferAnchors},
#' Default: 'LogNormalize'
#' @param reduction reduction method used during \code{FindTransferAnchors}, Default: 'cca'
#' @param method Method to use for annotation, either "Seurat" or "nmf".
#' @param min_prediction_score Minimum prediction score for the annotation to be considered valid.
#' Labels with a prediction score below this value will be labelled as "Unknown". Default: 0,
#' meaning that not threshold is used.
#' @param skip_normalization If TRUE, the datasets will not be normalized prior to annotation.
#' This parameter only has an effect if `method` is set to "nmf". The default layer will be
#' used, so use this setting at your own risk.
#' @param verbose If TRUE, print messages about the progress of the annotation.
#' @param ... Additional parameters to \code{FindTransferAnchors}
#'
#' @return The initial Seurat object with additional annotation columns added to its metadata.
#'
#' @details The "Seurat" method is a wrapper for the \code{FindTransferAnchors} and \code{TransferData}
#' functions from Seurat, followed by an optional summary per cluster.
#'
#' @examples
#' \dontrun{
#' # Download reference file
#' reference <-
#'   readRDS(url(
#'     paste0(
#'       "https://pixelgen-technologies-datasets.s3.eu-north-1.amazonaws.com/",
#'       "mpx-analysis/next/R/pbmc_annotation.rds"
#'     )
#'   ))
#'
#' # Does not work on the test data due to small size - ignore for now.
#' seur <- ReadPNA_Seurat(minimal_pna_pxl_file())
#'
#' seur <- AnnotateCells(seur,
#'   reference = reference, query_assay = "PNA"
#' )
#' }
#'
#' @rdname AnnotateCells
#'
#' @export
#'
AnnotateCells <- function(
  object,
  reference = reference,
  summarize_by_column = NULL,
  reference_assay = "ADT",
  query_assay = "PNA",
  reference_groups = c("celltype.l1", "celltype.l2"),
  normalization_method = c("LogNormalize", "SCT", "CLR"),
  reduction = "cca",
  method = c("Seurat", "nmf"),
  min_prediction_score = 0,
  skip_normalization = FALSE,
  verbose = TRUE,
  ...
) {
  expect_Seurat()

  # Input validation
  assert_class(object, "Seurat")
  assert_class(reference, "Seurat")
  assert_col_in_data(summarize_by_column, object[[]], allow_null = TRUE)
  assert_single_value(reference_assay, "string")
  assert_single_value(query_assay, "string")
  assert_vector(reference_groups, "character", n = 1)
  for (ref_group in reference_groups) {
    assert_col_in_data(ref_group, reference[[]])
  }
  assert_single_value(reduction, "string")
  method <- match.arg(method, c("Seurat", "nmf"))
  assert_single_value(skip_normalization, "bool")
  normalization_method <- match.arg(normalization_method, c("LogNormalize", "SCT", "CLR"))
  if (normalization_method == "CLR" && method == "Seurat") {
    cli::cli_abort(
      c("x" = "CLR normalization is not supported with the {.val Seurat} method.")
    )
  }
  assert_single_value(min_prediction_score, "numeric")
  assert_within_limits(min_prediction_score, c(0, 1))
  assert_single_value(verbose, "bool")

  groups <- reference[[]][, reference_groups, drop = FALSE] %>% as.list()

  object <- switch(method,
    "Seurat" = .seurat_annotaton(
      object,
      reference,
      query_assay,
      reference_assay,
      groups,
      normalization_method,
      skip_normalization,
      reduction,
      min_prediction_score,
      verbose,
      ...
    ),
    "nmf" = .seeded_nmf_annotaton(
      object,
      reference,
      query_assay,
      reference_assay,
      groups,
      normalization_method,
      skip_normalization,
      min_prediction_score,
      verbose
    )
  )

  if (!is.null(summarize_by_column)) {
    # Summarize the transferred cell types by retaining the most common cell type
    # per factor level in the user-provided column
    summarized_results <- lapply(reference_groups, function(ref_group) {
      object[[]] %>%
        group_by(!!sym(ref_group), !!sym(summarize_by_column)) %>%
        count() %>%
        group_by(!!sym(summarize_by_column)) %>%
        arrange(desc(n)) %>%
        slice_head(n = 1) %>%
        select(-n)
    }) %>%
      set_names(reference_groups)

    updated_metadata <- object[[]] %>%
      as_tibble(rownames = "component") %>%
      select(component, !!sym(summarize_by_column))

    for (ref_group in reference_groups) {
      updated_metadata <- updated_metadata %>%
        left_join(summarized_results[[ref_group]], by = summarize_by_column) %>%
        rename_with(~ paste0(.x, "_summary"), !!sym(ref_group))
    }

    updated_metadata <- updated_metadata %>%
      data.frame(row.names = "component") %>%
      select(-!!sym(summarize_by_column))

    object <- SeuratObject::AddMetaData(
      object = object,
      metadata = updated_metadata
    )
  }

  return(object)
}

#' Run cell type annotation using Seurat method
#'
#' @param object Seurat object containing the data to be annotated.
#' @param reference Seurat object containing the reference data for annotation.
#' @param query_assay Name of the assay in the query object to be annotated.
#' @param reference Assay in the reference object to be used for annotation.
#' @param groups List of character vectors specifying the groupings in the reference object.
#' @param normalization_method Method for normalizing the data, e.g., "LogNormalize".
#' @param skip_normalization Logical indicating whether to skip normalization.
#' @param reduction Type of dimensionality reduction to use, e.g., "cca".
#' @param min_prediction_score Minimum prediction score for cell type assignment.
#' @param verbose Logical indicating whether to print progress messages.
#' @param ... Additional arguments passed to `Seurat::FindTransferAnchors`.
#'
#' @noRd
#'
.seurat_annotaton <- function(
  object,
  reference,
  query_assay,
  reference_assay,
  groups,
  normalization_method,
  skip_normalization,
  reduction,
  min_prediction_score,
  verbose,
  ...
) {
  # change assay classes if needed
  reference_class <- class(reference[[reference_assay]])
  query_class <- class(object[[query_assay]])
  assay_temp <- object[[query_assay]]
  if (reference_class != query_class) {
    if (query_class %in% c("PNAAssay", "PNAAssay5")) {
      new_assay_class <- switch(query_class,
        "PNAAssay" = "Assay",
        "PNAAssay5" = "Assay5"
      )
      assay_temp <- try(
        {
          as(assay_temp, new_assay_class)
        },
        silent = TRUE
      )
      if (inherits(assay_temp, "try-error")) {
        cli::cli_abort(
          c("x" = "Failed to convert query assay class {.cls {query_class}} to {.cls new_assay_class}")
        )
      }
    }
    if (class(assay_temp) != reference_class) {
      reference[[reference_assay]] <- switch(new_assay_class,
        "Assay" = SeuratObject::CreateAssayObject(counts = LayerData(reference[[reference_assay]], layer = "counts")),
        "Assay5" = SeuratObject::CreateAssay5Object(counts = LayerData(reference[[reference_assay]], layer = "counts"))
      )
    }
    if (class(assay_temp) != class(reference[[reference_assay]])) {
      cli::cli_abort(
        c(
          "i" = paste0(
            "Cannot combine reference assay class {.cls {reference_class}} ",
            "with query assay class {.cls {query_class}}. "
          ),
          "x" = "Assay class mismatch."
        )
      )
    }
  }
  object_temp <- object
  object_temp[[query_assay]] <- assay_temp

  # Set variable features if missing
  if (is.null(VariableFeatures(reference))) {
    VariableFeatures(reference) <- rownames(reference[[reference_assay]])
  }
  if (is.null(VariableFeatures(object_temp))) {
    VariableFeatures(object_temp) <- rownames(assay_temp)
  }

  # Apply normalization
  if (!skip_normalization) {
    if (verbose && check_global_verbosity()) {
      cli::cli_alert_info("Normalizing reference and query data using method {.str {normalization_method}}")
    }
    reference <-
      Seurat::NormalizeData(
        reference,
        assay = reference_assay,
        normalization.method = normalization_method,
        margin = 2,
        verbose = verbose
      )
    object_temp <-
      Seurat::NormalizeData(
        object_temp,
        assay = query_assay,
        normalization.method = normalization_method,
        margin = 2,
        verbose = verbose
      )
  }

  # Find the anchors to bridge the two datasets
  anchors <- Seurat::FindTransferAnchors(
    reference = reference,
    query = object_temp,
    reference.assay = reference_assay,
    query.assay = query_assay,
    normalization.method = normalization_method,
    reduction = reduction,
    l2.norm = TRUE,
    verbose = verbose,
    ...
  )

  # Transfer the data and keep the celltypes
  predictions <- Seurat::TransferData(
    anchorset = anchors,
    refdata = groups,
    weight.reduction = reduction
  )

  # Format results
  annotations <- lapply(predictions, function(x) x[, "predicted.id"]) %>% as.data.frame()

  if (min_prediction_score > 0) {
    prediction_scores <- lapply(predictions, function(x) x[, "prediction.score.max"]) %>% as.data.frame()
    annotations[prediction_scores < min_prediction_score] <- "Unknown"
  }

  if (all(names(annotations) == c("celltype.l1", "celltype.l2"))) {
    colnames(annotations) <- c("l1_annotation", "l2_annotation")
  } else {
    colnames(annotations) <- paste0(colnames(annotations), "_annotation")
  }

  # Add the transferred celltypes to the object metadata
  object <- SeuratObject::AddMetaData(
    object = object,
    metadata = annotations
  )

  return(object)
}


#' Run cell type annotation using seeded NMF
#'
#' @param object Seurat object containing the data to be annotated.
#' @param reference Seurat object containing the reference data for annotation.
#' @param query_assay Name of the assay in the query object to be annotated.
#' @param reference_assay Assay in the reference object to be used for annotation.
#' @param groups List of character vectors specifying the groupings in the reference object.
#' @param normalization_method Method for normalizing the data, e.g., "LogNormalize".
#' @param skip_normalization Logical indicating whether to skip normalization of the data.
#' @param min_prediction_score Minimum prediction score for cell type assignment.
#' @param verbose Logical indicating whether to print progress messages.
#'
#' @noRd
#'
.seeded_nmf_annotaton <- function(
  object,
  reference,
  query_assay,
  reference_assay,
  groups,
  normalization_method,
  skip_normalization,
  min_prediction_score,
  verbose
) {
  expect_RcppML()

  # Apply normalization
  if (!skip_normalization) {
    if (verbose && check_global_verbosity()) {
      cli::cli_alert_info("Normalizing reference and query data using method {.str {normalization_method}}")
    }
    reference <-
      Seurat::NormalizeData(
        reference,
        normalization.method = normalization_method,
        margin = 2,
        verbose = verbose
      )
    object <-
      Seurat::NormalizeData(
        object,
        normalization.method = normalization_method,
        margin = 2,
        verbose = verbose
      )
  }

  if (verbose && check_global_verbosity()) {
    cli::cli_alert_info("Fetching overlapping data from reference and query assays.")
  }

  # Get data
  reference_data <- SeuratObject::LayerData(reference, assay = reference_assay)
  target_data <- SeuratObject::LayerData(object, assay = query_assay)

  # Subset data to use overlapping markers
  shared_markers <- intersect(rownames(reference_data), rownames(target_data))
  reference_data <- reference_data[shared_markers, ]
  target_data <- target_data[shared_markers, ]

  if (verbose && check_global_verbosity()) {
    cli::cli_alert_info(
      c(
        "Running seeded NMF on query cells for cell type annotations {.str {names(groups)}} ",
        "in reference data."
      )
    )
  }

  predicted_annotations <- lapply(names(groups), function(nm) {
    W <- .get_w_matrix(
      target_data = target_data,
      reference_data = reference_data,
      nCells_per_group = 100,
      minCells_per_celltype = 10,
      groups = groups[[nm]],
      seed = 123,
      verbose = FALSE
    )

    # Run NNLS
    proj_expr <- RcppML::project(as(target_data, "dgCMatrix"), W, L1 = 0.01)

    # Convert predicted values to proportions
    prop <- apply(proj_expr, 2, function(x) {
      prop.table(x)
    })
    rownames(prop) <- colnames(W)
    colnames(prop) <- colnames(target_data)

    # Predicted cell type
    predicted_celltype <- apply(prop, 2, function(x) {
      rownames(prop)[which.max(x)]
    })

    if (min_prediction_score > 0) {
      prediction_score <- apply(prop, 2, max)
      predicted_celltype[prediction_score < min_prediction_score] <- "Unknown"
    }

    return(predicted_celltype)
  }) %>% set_names(nm = names(groups))

  if (verbose && check_global_verbosity()) {
    cli::cli_alert_info("Labelling cells.")
  }

  # Format results
  predicted_annotations <- predicted_annotations %>% as.data.frame()

  if (all(names(predicted_annotations) == c("celltype.l1", "celltype.l2"))) {
    colnames(predicted_annotations) <- c("l1_annotation", "l2_annotation")
  } else {
    colnames(predicted_annotations) <- paste0(colnames(predicted_annotations), "_annotation")
  }

  # Add the transferred celltypes to the object metadata
  object <- SeuratObject::AddMetaData(
    object = object,
    metadata = predicted_annotations
  )

  if (verbose && check_global_verbosity()) {
    cli::cli_alert_success("Finished")
  }

  return(object)
}

#' Compute an enrichment score for each marker and cell type
#'
#' @param x A matrix of protein abundance values, where rows are proteins
#' and columns are cells.
#' @param groups A vector indicating the group of each column in `x`.
#'
#' @return A matrix where rows are proteins and columns are groups, with
#' enrichment scores calculated as the ratio of the mean expression
#' of each protein in a group to the mean expression of that protein
#'
#' @noRd
#'
.get_abundance_profiles <- function(
  x,
  groups
) {
  # Calculate means
  row_means <- do.call(cbind, lapply(unique(groups), function(grp) {
    Matrix::rowMeans(x[, groups == grp, drop = FALSE])
  }))
  colnames(row_means) <- unique(groups)

  # Calculate enrichment scores
  W <- do.call(bind_cols, lapply(colnames(row_means), function(grp) {
    x1 <- row_means[, grp]
    x2 <- row_means[, -which(colnames(row_means) == grp)]
    x2 <- Matrix::rowMeans(x2) + 1
    y <- tibble(x1 / x2) %>% set_names(nm = grp)
    return(y)
  })) %>%
    as.matrix()
  rownames(W) <- rownames(x)

  W <- apply(W, 2, function(x) {
    x / max(x)
  })

  return(W)
}

#' Compute an "enrichment" matrix to be used as a seed for NMF
#'
#' @param target_data A matrix of protein abundance values for the target dataset.
#' @param reference_data A matrix of protein abundance values for the reference dataset.
#' @param nCells_per_group Maximum number of cells to sample per group in the `reference_data`.
#' @param minCells_per_celltype Minimum number of cells required per cell type in the
#' `reference_data`.
#' @param groups A vector group labels for each each column in the `reference_data`.
#' @param seed Random seed for reproducibility.
#' @param verbose Logical indicating whether to print progress messages.
#'
#' @return A matrix with enrichment scores where rows are proteins and columns are groups.
#'
#' @noRd
#'
.get_w_matrix <- function(
  target_data,
  reference_data,
  nCells_per_group = 50,
  minCells_per_celltype = 10,
  groups,
  seed = 123,
  verbose = TRUE
) {
  # Rescale data to ensure positive values
  if (any(target_data < 0)) {
    if (verbose) warn("Found negative values in input matrix")
    if (verbose) cli_alert_info("Rescaling data to ensure positive values")
    object <- t(apply(object, 1, function(x) (x - min(x)) / (max(x) - min(x))))
  }

  # Sample cells
  if (verbose) {
    cli_alert(
      glue(
        "  Downsampling data to include a maximum of ",
        "{nCells_per_group} cells per cell type"
      )
    )
  }

  set.seed(seed)
  components <- tibble(component = colnames(reference_data), group = groups) %>%
    group_by(group) %>%
    slice(sample(min(nCells_per_group, n())))

  # Check that all celltypes have enough cells
  n_cells <- components %>% count()
  if (any(n_cells$n < minCells_per_celltype)) {
    too_small <- n_cells %>%
      filter(n < minCells_per_celltype)
    if (verbose) {
      cli_alert(
        c(
          "  Cell type{?s} {.val {too_small$group}} ",
          "have too few cells and will be excluded"
        )
      )
    }
    components <- components %>%
      filter(!group %in% too_small$group)
  }
  if (verbose) cli_alert("  Kept {.val {length(unique(components$group))}} cell types after filtering")
  components <- components %>% pull(component)
  new_groups <- set_names(groups, nm = colnames(reference_data))[components]

  # Compute cell type abundance profiles
  if (verbose) cli_alert("  Calculating cell type abundance profiles")
  W <- .get_abundance_profiles(x = reference_data[, components], groups = new_groups)

  return(W)
}
