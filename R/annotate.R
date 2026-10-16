#' Automatic annotation of cell types
#'
#' This function transfers cell type annotations from a reference Seurat dataset
#' to a query Seurat dataset using non-negative least squares (\code{"nnls"}, default)
#' or Seurat anchor transfer (\code{"Seurat"}). For each column specified in
#' \code{reference_groups}, a corresponding column of predicted cell type annotations
#' and a \code{<group>_score} column of per-cell prediction scores are added to the
#' metadata of the query dataset. In addition, one can optionally
#' summarize these annotations for a specific clustering solution, which adds a
#' \code{<group>_summary} column containing the most common annotation per cluster.
#'
#' @param object Seurat object to which you want to add celltype annotation
#' @param reference Seurat object that contains annotated reference data, see
#' \code{\link{read_pbmc_reference}} and examples for obtaining a PBMC reference dataset
#' @param summarize_by_column If NULL, return only the cell type annotation per cell.
#' If the name of a column in the object metadata is provided, the annotations will
#' be summarized to the factors of this column, i.e. the most common annotation per
#' factor level will be retained. Useful for annotating the output of the Seurat
#' clustering, for example, Default: NULL
#' @param reference_assay Assay in the reference dataset used for annotation,
#' Default: 'PNA'
#' @param query_assay Assay in the query dataset used for annotation,
#' Default: 'PNA'
#' @param reference_groups A character vector of reference metadata columns to use for
#' annotation, Default: \code{c("celltype.l1", "celltype.l2")}
#' @param normalization_method Normalization method used prior to annotation (with margin = 2),
#' either "CLR", "LogNormalize", or "SCT", Default: 'CLR'
#' @param reduction Reduction method used during \code{FindTransferAnchors} when
#' \code{method = "Seurat"}, Default: 'cca'
#' @param method Method to use for annotation, either "nnls", "Seurat", or "nmf"
#' (alias for "nnls"), Default: 'nnls'
#' @param min_prediction_score Minimum prediction score for the annotation to be considered valid.
#' Labels with a prediction score below this value will be labelled as "Unknown", Default: 0,
#' meaning that no threshold is used.
#' @param skip_normalization If TRUE, the datasets will not be normalized prior to annotation.
#' The default layer will be used, so use this setting at your own risk, Default: FALSE
#' @param verbose If TRUE, print messages about the progress of the annotation, Default: TRUE
#' @param ... Additional parameters passed to \code{FindTransferAnchors} when \code{method = "Seurat"}
#'
#' @return The initial Seurat object with additional annotation columns added to its metadata.
#'
#' @details The "nnls" (or "nmf") method projects query cells onto feature-weighted reference
#' cell type profiles using non-negative least squares (\code{RcppML::project}), followed by
#' local fine-tuning across the top candidate cell types. The "Seurat" method is a wrapper for
#' the \code{FindTransferAnchors} and \code{TransferData} functions from Seurat. Both methods
#' support an optional summary per cluster.
#'
#' Merged \pkg{SeuratObject} v5 assays (\code{Assay5}, including \code{PNAAssay5}) must not
#' retain split merge layers (\code{counts.1}, \code{data.1}, etc.).
#' Run \code{JoinLayers()} on the \code{Seurat} object (or on the assay) after \code{merge()}
#' before calling \code{AnnotateCells()}.
#'
#' @examples
#' \dontrun{
#' reference <- read_pbmc_reference()
#'
#' seur <- ReadPNA_Seurat(minimal_pna_pxl_file())
#'
#' seur <- AnnotateCells(
#'   seur,
#'   reference = reference
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
  reference_assay = "PNA",
  query_assay = "PNA",
  reference_groups = c("celltype.l1", "celltype.l2"),
  normalization_method = c("CLR", "LogNormalize", "SCT"),
  reduction = "cca",
  method = c("nnls", "Seurat", "nmf"),
  min_prediction_score = 0,
  skip_normalization = FALSE,
  verbose = TRUE,
  ...
) {
  expect_Seurat()

  # Input validation
  assert_class(object, "Seurat")
  assert_class(reference, "Seurat")

  .assert_assay_layers_joined(object[[query_assay]], "query", query_assay)
  .assert_assay_layers_joined(reference[[reference_assay]], "reference", reference_assay)

  assert_col_in_data(summarize_by_column, object[[]], allow_null = TRUE)
  assert_single_value(reference_assay, "string")
  assert_single_value(query_assay, "string")
  assert_vector(reference_groups, "character", n = 1)
  for (ref_group in reference_groups) {
    assert_col_in_data(ref_group, reference[[]])
  }
  assert_single_value(reduction, "string")
  method <- match.arg(method, c("nnls", "Seurat", "nmf"))
  if (method == "nmf") {
    method <- "nnls"
  }
  assert_single_value(skip_normalization, "bool")
  normalization_method <- match.arg(normalization_method, c("CLR", "LogNormalize", "SCT"))
  if (normalization_method == "CLR" && method == "Seurat" && !skip_normalization) {
    cli::cli_abort(
      c("x" = "CLR normalization is not supported with the {.val Seurat} method when {.arg skip_normalization} is {.val FALSE}.")
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
    "nnls" = .seeded_nmf_annotaton(
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
  # Seurat::FindTransferAnchors only accepts "LogNormalize" or "SCT" for normalization.method
  anchor_norm_method <- if (normalization_method == "CLR") "LogNormalize" else normalization_method
  anchors <- Seurat::FindTransferAnchors(
    reference = reference,
    query = object_temp,
    reference.assay = reference_assay,
    query.assay = query_assay,
    normalization.method = anchor_norm_method,
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
  if (!is.data.frame(predictions)) {
    annotations <- lapply(predictions, function(x) x[, "predicted.id"]) %>% data.frame(check.names = FALSE)
    prediction_scores <- lapply(predictions, function(x) x[, "prediction.score.max"]) %>% data.frame(check.names = FALSE)
  } else {
    annotations <- data.frame(predictions[, "predicted.id"], row.names = rownames(predictions), check.names = FALSE)
    prediction_scores <- data.frame(predictions[, "prediction.score.max"], row.names = rownames(predictions), check.names = FALSE)
    colnames(annotations) <- names(groups)
    colnames(prediction_scores) <- names(groups)
  }

  if (min_prediction_score > 0) {
    annotations[prediction_scores < min_prediction_score] <- "Unknown"
  }

  colnames(prediction_scores) <- paste0(names(groups), "_score")
  metadata_to_add <- cbind(annotations, prediction_scores)

  # Add the transferred celltypes and prediction scores to the object metadata
  object <- SeuratObject::AddMetaData(
    object = object,
    metadata = metadata_to_add
  )

  return(object)
}


#' Fail fast when an Assay5 still has split merge layers (counts.1, data.1, ...).
#'
#' @param assay Assay from `object[[assay_name]]`.
#' @param role `"query"` or `"reference"` for the error message.
#' @param assay_name Assay name for the error message.
#'
#' @noRd
#'
.assert_assay_layers_joined <- function(assay, role = " ", assay_name) {
  if (!methods::is(assay, "Assay5")) {
    return(invisible(NULL))
  }
  lys <- SeuratObject::Layers(assay)
  split_layers <- grep("^(counts|data)\\.\\d+$", lys, value = TRUE, perl = TRUE)
  if (length(split_layers) == 0L) {
    return(invisible(NULL))
  }
  cli::cli_abort(
    c(
      "x" = "The {role}assay {.val {assay_name}} has unjoined v5 layers ({.val {split_layers}}).",
      "i" = "Call {.fn JoinLayers} on the {.cls Seurat} object (or on that assay)
      after {.fn merge}, then retry {.fn AnnotateCells} with {.arg method} = {.val nmf}."
    )
  )
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

  # When normalization is required for NMF: normalize the reference in place, and build
  # a separate minimal Seurat from the query counts for NormalizeData.
  if (!skip_normalization) {
    if (verbose && check_global_verbosity()) {
      cli::cli_alert_info("Normalizing reference and query data using method {.str {normalization_method}}")
    }
    reference <- Seurat::NormalizeData(
      reference,
      assay = reference_assay,
      normalization.method = normalization_method,
      margin = 2,
      verbose = verbose
    )
    query_counts <- SeuratObject::LayerData(object, assay = query_assay, layer = "counts")
    query_for_nm <- SeuratObject::CreateSeuratObject(
      counts = query_counts,
      assay = query_assay,
      meta.data = data.frame(row.names = colnames(query_counts))
    )
    query_for_nm <- Seurat::NormalizeData(
      query_for_nm,
      assay = query_assay,
      normalization.method = normalization_method,
      margin = 2,
      verbose = verbose
    )
  }

  if (verbose && check_global_verbosity()) {
    cli::cli_alert_info("Fetching overlapping data from reference and query assays.")
  }

  reference_data <- SeuratObject::LayerData(reference, assay = reference_assay)
  target_data <- if (!skip_normalization) {
    SeuratObject::LayerData(query_for_nm, assay = query_assay)
  } else {
    SeuratObject::LayerData(object, assay = query_assay)
  }

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
    w_res <- .get_w_matrix(
      reference_data = reference_data,
      nCells_per_group = 500,
      minCells_per_celltype = 5,
      groups = groups[[nm]],
      seed = 123,
      verbose = FALSE
    )

    W <- w_res$W
    rm <- w_res$row_means
    w_var <- w_res$within_var
    fw <- w_res$feature_weights
    parent_map <- w_res$parent_map
    u_parents <- unique(parent_map)

    # Scale query features by global root-F weights and run Stage 1 NNLS
    weighted_target <- as(sweep(target_data, 1, fw, "*"), "dgCMatrix")
    proj_expr <- RcppML::project(weighted_target, W, L1 = 0.05)

    # Cells with a very low total projection across all cell types are labeled
    # as "Unknown" to avoid spurious annotations.
    failed_cells <- colSums(proj_expr) < 1e-2

    # Add a small value to the scores of failed cells to avoid issues with zero values in downstream calculations.
    proj_expr[, failed_cells] <- proj_expr[, failed_cells] + 1e-8

    # Aggregate Stage 1 projections by parent group to compute prediction proportions
    parent_proj <- vapply(u_parents, function(p) {
      sub_idx <- which(parent_map == p)
      if (length(sub_idx) == 1L) {
        proj_expr[sub_idx, ]
      } else {
        colSums(proj_expr[sub_idx, , drop = FALSE])
      }
    }, numeric(ncol(proj_expr)))
    parent_proj <- t(parent_proj)
    rownames(parent_proj) <- u_parents

    # Convert Stage 1 predicted values to proportions
    prop <- apply(parent_proj, 2, function(x) {
      prop.table(x)
    })
    if (is.null(dim(prop))) {
      prop <- matrix(prop, nrow = 1L, dimnames = list(u_parents, colnames(target_data)))
    } else {
      rownames(prop) <- u_parents
    }

    # Stage 2: Pairwise local fine-tuning on top-2 candidate parent groups per cell
    # (each represented by its highest-scoring sub-centroid for that cell, and only
    # when the runner-up has non-trivial Stage 1 support relative to the top candidate)
    k_cands <- min(2L, length(u_parents))
    min_rel_score <- 0.25
    sub_indices_by_parent <- split(seq_along(parent_map), parent_map)[u_parents]
    predicted_celltype <- vapply(seq_len(ncol(target_data)), function(i) {
      if (k_cands <= 1L) {
        return(u_parents[1L])
      }
      scores_i <- proj_expr[, i]
      best_sub_idx <- vapply(sub_indices_by_parent, function(idxs) {
        if (length(idxs) == 1L) idxs else idxs[which.max(scores_i[idxs])]
      }, integer(1))
      sub_scores <- scores_i[best_sub_idx]
      ord <- order(sub_scores, decreasing = TRUE)
      top_k <- ord[seq_len(k_cands)]
      top_k <- top_k[sub_scores[top_k] >= sub_scores[top_k[1L]] * min_rel_score]
      if (length(top_k) <= 1L) {
        return(u_parents[ord[1L]])
      }
      cands <- colnames(W)[best_sub_idx[top_k]]
      sub_rm <- rm[, cands, drop = FALSE]

      # Weight features by local between-candidate variance relative to pooled within-group variance,
      # capped at the 90th percentile so a single batch-shifted marker cannot override the consensus
      local_b_var <- apply(sub_rm, 1, stats::var)
      local_fw <- sqrt(local_b_var / w_var)
      if (all(local_fw == 0)) {
        return(unname(parent_map[cands[1L]]))
      }
      cap_val <- stats::quantile(local_fw[local_fw > 0], 0.90, names = FALSE)
      local_fw <- pmin(local_fw, cap_val)

      sub_W <- sweep(sub_rm, 1, local_fw, "*")
      sub_W <- sweep(sub_W, 2, sqrt(colSums(sub_W^2)) + 1e-8, "/")
      sub_q <- local_fw * target_data[, i]
      best_cand <- cands[which.max(crossprod(sub_W, sub_q))]
      unname(parent_map[best_cand])
    }, character(1))

    prediction_score <- apply(prop, 2, max)
    prediction_score[failed_cells] <- 0
    names(prediction_score) <- colnames(target_data)

    if (min_prediction_score > 0) {
      predicted_celltype[prediction_score < min_prediction_score] <- "Unknown"
    }

    predicted_celltype[failed_cells] <- "Unknown"
    names(predicted_celltype) <- colnames(target_data)

    out <- data.frame(
      celltype = predicted_celltype,
      score = prediction_score,
      row.names = colnames(target_data),
      check.names = FALSE
    )
    colnames(out) <- c(nm, paste0(nm, "_score"))
    return(out)
  })

  if (verbose && check_global_verbosity()) {
    cli::cli_alert_info("Labelling cells.")
  }

  # Format results
  predicted_annotations <- do.call(cbind, predicted_annotations)

  # Add the transferred celltypes and prediction scores to the object metadata
  object <- SeuratObject::AddMetaData(
    object = object,
    metadata = predicted_annotations
  )

  if (verbose && check_global_verbosity()) {
    cli::cli_alert_success("Finished")
  }

  return(object)
}

#' Compute feature-weighted reference profiles for NNLS annotation
#'
#' @param x A matrix of protein abundance values, where rows are proteins
#' and columns are cells.
#' @param groups A vector indicating the group of each column in `x`.
#'
#' @return A list containing:
#' \itemize{
#'   \item `W`: L2-normalized, feature-weighted centroid matrix (proteins x groups)
#'   \item `row_means`: Unweighted group centroid matrix (proteins x groups)
#'   \item `within_var`: Pooled within-group variance vector per protein
#'   \item `feature_weights`: Root-F feature weight vector per protein
#' }
#'
#' @noRd
#'
.get_abundance_profiles <- function(
  x,
  groups,
  max_subclusters = 3L,
  min_subcluster_size = 12L
) {
  groups <- as.character(groups)
  ugroups <- unique(groups)

  # Split heterogeneous/large groups into up to `max_subclusters` sub-centroids
  # so broad lineages (e.g. L1 Mono = Classical + Non-classical Mono) are represented
  # by their constituent modes rather than a single blurred average
  sub_groups <- character(length(groups))
  parent_map <- character(0)

  for (grp in ugroups) {
    idx <- which(groups == grp)
    n_grp <- length(idx)
    if (n_grp >= min_subcluster_size * 2L && max_subclusters > 1L) {
      k <- min(max_subclusters, floor(n_grp / min_subcluster_size))
      sub_x <- as.matrix(x[, idx, drop = FALSE])
      # Remove zero-variance rows before k-means
       cl <- tryCatch({
        km <- stats::kmeans(t(sub_x), centers = k, nstart = 5L, iter.max = 20L)
        km$cluster
      }, error = function(e) {
        rep(1L, n_grp)
      })
      cl_tab <- table(cl)
      largest_cl <- as.integer(names(cl_tab)[which.max(cl_tab)])
      for (small_c in as.integer(names(cl_tab)[cl_tab < 5L])) {
        cl[cl == small_c] <- largest_cl
      }
      if (length(unique(cl)) > 1L) {
        sub_ids <- paste0(grp, "__sub", cl)
        sub_groups[idx] <- sub_ids
        for (sid in unique(sub_ids)) {
          parent_map[sid] <- grp
        }
      } else {
        sub_groups[idx] <- grp
        parent_map[grp] <- grp
      }
    } else {
      sub_groups[idx] <- grp
      parent_map[grp] <- grp
    }
  }

  usub <- unique(sub_groups)

  # Calculate robust sub-group centroids using a 15% trimmed mean across cells
  row_means <- vapply(usub, function(grp) {
    sub_mat <- as.matrix(x[, sub_groups == grp, drop = FALSE])
    apply(sub_mat, 1, base::mean, trim = 0.15)
  }, numeric(nrow(x)))
  colnames(row_means) <- usub
  rownames(row_means) <- rownames(x)

  # Pooled within-sub-group variance per marker (computed directly on sparse/dense columns)
  within_var <- vapply(usub, function(grp) {
    sub_x <- x[, sub_groups == grp, drop = FALSE]
    n_grp <- ncol(sub_x)
    if (n_grp <= 1L) {
      return(rep(0, nrow(sub_x)))
    }
    grp_mean <- Matrix::rowMeans(sub_x)
    pmax(0, (Matrix::rowSums(sub_x^2) - n_grp * grp_mean^2) / (n_grp - 1L))
  }, numeric(nrow(x)))
  within_var <- rowMeans(within_var, na.rm = TRUE) + 0.05

  # Between-group variance per marker and root-F feature weights (capped at 95th percentile
  # for robustness against single-marker titration shifts across batches)
  between_var <- if (ncol(row_means) > 1L) {
    apply(row_means, 1, stats::var)
  } else {
    rep(1, nrow(row_means))
  }
  feature_weights <- sqrt(between_var / within_var)
  if (all(feature_weights == 0)) {
    feature_weights <- rep(1, nrow(row_means))
  } else {
    cap_fw <- stats::quantile(feature_weights[feature_weights > 0], 0.95, names = FALSE)
    feature_weights <- pmin(feature_weights, cap_fw)
  }

  # Weight group centroids by feature weights and L2-normalize columns
  W <- sweep(row_means, 1, feature_weights, "*")
  W <- sweep(W, 2, sqrt(colSums(W^2)) + 1e-8, "/")

  return(list(
    W = W,
    row_means = row_means,
    within_var = within_var,
    feature_weights = feature_weights,
    parent_map = parent_map
  ))
}

#' Compute feature-weighted reference matrices to be used for seeded NMF / NNLS
#'
#' @param reference_data A matrix of protein abundance values for the reference dataset.
#' @param nCells_per_group Maximum number of cells to sample per group in the `reference_data`.
#' @param minCells_per_celltype Minimum number of cells required per cell type in the
#' `reference_data`.
#' @param groups A vector group labels for each column in the `reference_data`.
#' @param seed Random seed for reproducibility.
#' @param verbose Logical indicating whether to print progress messages.
#'
#' @return A list returned by `.get_abundance_profiles()` containing `W`,
#' `row_means`, `within_var`, and `feature_weights`.
#'
#' @noRd
#'
.get_w_matrix <- function(
  reference_data,
  nCells_per_group = 500,
  minCells_per_celltype = 5,
  groups,
  seed = 123,
  verbose = TRUE
) {
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
  kept_components <- components %>% pull(component)
  new_groups <- set_names(groups, nm = colnames(reference_data))[kept_components]

  # Compute cell type abundance profiles
  if (verbose) cli_alert("  Calculating cell type abundance profiles")
  w_res <- .get_abundance_profiles(x = reference_data[, kept_components, drop = FALSE], groups = new_groups)

  return(w_res)
}


#' Read PBMC cell annotation reference
#'
#' This function reads the cell annotation reference data from the pixelatorR package's extdata directory. This data
#' contains 720 annotated PBMC components with their corresponding counts of 158 proteins. The annotation has been
#' performed at two levels of granularity: "celltype.l1" (broad categories) and "celltype.l2" (more specific subtypes).
#' When loading the data through this function, a Seurat object is created with the counts data stored in the "PNA"
#' assay and the cell annotations included in the metadata.
#'
#' @return A Seurat object containing the cell annotation and counts data.
#'
#' @export
#'
read_pbmc_reference <-
  function() {
    cell_annotation <-
      system.file(
        "extdata", "PBMC_cell_annotation.csv",
        package = "pixelatorR"
      ) %>%
      read.csv() %>%
      column_to_rownames("sample_component")

    counts <-
      system.file(
        "extdata", "PBMC_cell_counts.csv",
        package = "pixelatorR"
      ) %>%
      read.csv(check.names = FALSE) %>%
      column_to_rownames("sample_component")

    counts <-
      counts[rownames(cell_annotation), ] %>%
      t() %>%
      as("dgCMatrix")

    # Create Seurat object
    CreateSeuratObject(
      counts = counts,
      assay = "PNA",
      meta.data = cell_annotation
    )
  }
