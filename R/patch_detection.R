#' Supervised patch detection
#'
#' A patch is defined as a subgraph which is enriched for a set of patch-specific protein
#' markers. A patch should typically have a different origin than the bulk of the PNA graph.
#' A typical example of a patch is a small piece of another cell, e.g. a patch of a B cell on
#' a T cell (receiver).
#'
#' This analysis tool requires a predefined set of patch-specific protein markers and is therefore
#' a supervised method. Patch-specific means that these markers are high abundant on the patch
#' and low abundant on the receiver cell. Optimal patch markers are those that are high abundant
#' and highly specific. The method is sensitive to the choice of patch markers and therefore
#' requires careful selection. See \code{\link{identify_markers_for_patch_analysis}} for more
#' information on how to select patch markers.
#'
#' Patch detection can be used to find multiple patches in a single PNA graph. These subgraphs
#' can be leveraged to study the patch protein composition, patch proximity scores, number of patches
#' and the fraction of the receiver cell covered by patches.
#' Note that patches can appear from cell debris or from artificial bleedover and doesn't
#' necessarily indicate a cell-cell interaction. If you aim to study patch composition in response to
#' a specific treatment or condition, it is recommended to include a control population for reference.
#'
#' Two algorithms are available for patch detection, described in the sections below.
#'
#' @section Expand and contract:
#' 1. Initialization: Start with the set of nodes `P` labelled by patch-specific protein markers.
#' 2. Expansion:
#'    - Expand `P` to include nodes that connect with at least 2 nodes in `P` within
#'      a 2-step neighborhood.
#'    - Expand `P` to include nodes that connect with at least 2 nodes in `P` within
#'      a 1-step neighborhood.
#' 3. Contract `P` by only keeping nodes with the highest in- out-degree ratio. The number of
#' kept nodes is controlled by the `contraction` parameter.
#' 4. (optional) Remove nodes from `P` at the patch border with low patch connectivity
#' 5. Construct the patch graph from `P` and split it into its connected components and remove
#' components smaller than `patch_nodes_threshold`
#' 6. (optional) Run an iterative community detection (Leiden) to split up weakly connected
#' patch components and repeat the `patch_nodes_threshold` filtering step. The resolution
#' is defined by `leiden_resolution` where a higher value is more likely to result in more
#' patches and vice versa.
#' 7. Label nodes in the original graph with the patch information. The largest
#' "patch" is labeled as 0 and should  correspond to the receiver cell graph. The rest
#' of the patches are labeled as 1, 2, etc. and correspond to patches ordered by
#' decreasing size.
#'
#' @section Local G:
#' 1. Run local G using the `patch_markers` UMI counts.
#' 2. Define `P` as the set of nodes with a p value below `pval_threshold`.
#' 3. (optional) Remove nodes from `P` at the patch border with low patch connectivity
#' 4. Construct the patch graph from `P` and split it into its connected components and remove
#' components smaller than `patch_nodes_threshold`.
#' 5. (optional) Run an iterative community detection (Leiden) to split up weakly connected
#' patch components and repeat the `patch_nodes_threshold` filtering step. The resolution
#' is defined by `leiden_resolution` where a higher value is more likely to result in more
#' patches and vice versa.
#' 6. Label nodes in the original graph with the patch information. The largest
#' "patch" is labeled as 0 and should  correspond to the receiver cell graph. The rest
#' of the patches are labeled as 1, 2, etc. and correspond to patches ordered by
#' decreasing size.
#'
#' @param cg A \code{CellGraph} object with PNA data.
#' @param patch_markers A character vector with the names of the markers that are
#' exclusively found on the patches, e.g. "CD41" for platelets.
#' @param receiver_markers An optional character vector with the names of the markers that are
#' exclusively found on the receiver cell.
#' @param k The number of nearest neighbors to consider for the expansion or for local G.
#' @param leiden_resolution The resolution parameter for the Leiden algorithm.
#' @param patch_nodes_threshold The minimum number of nodes to consider a patch.
#' @param prune_patch_edge A logical indicating if edge nodes of the patches should be pruned.
#' @param leiden_refinement A logical indicating if the patch graph should be refined into
#' smaller communities using Leiden. This is useful if there are weakly connected patches.
#' @param method A character indicating the method to use for patch detection. Must be one
#' of "expand_contract" or "local_g".
#' @param contraction A numeric value between 0 and 1 that controls the number of nodes to keep
#' in the "expand_contract" method. Higher values will increase contraction and result in smaller
#' patches.
#' @param pval_threshold A numeric value that controls the p value threshold for the
#' "local_g" method.
#' @param seed Set seed for reproducibility
#' @param verbose A logical indicating if messages should be printed to the console.
#'
#' @examples
#' library(tidygraph)
#' library(dplyr)
#' library(ggplot2)
#' library(Matrix)
#'
#' se <- ReadPNA_Seurat(minimal_pna_pxl_file()) %>%
#'   LoadCellGraphs(cells = colnames(.)[1], add_layouts = TRUE)
#'
#' cg <- CellGraphs(se)[[1]]
#'
#' protein_props <- cg@counts %>% Matrix::colSums() %>% prop.table()
#' patch_props <- protein_props
#' patch_props["CD8"] <- 0.4
#' patch_props <- patch_props %>% prop.table()
#'
#' # Sample node indices for a. new patch
#' inds <- 1
#' patch_size <- 1000
#' xyz <- cg@layout$wpmds_3d %>% as.matrix()
#' xyz_center <- xyz[inds, , drop = FALSE]
#' dists <- 1 - cos_dist(A = xyz_center, B = xyz) %>% as.vector()
#' inds_replace <- order(dists)[1:patch_size]
#' counts_with_patch <- cg@counts
#'
#' # Create a count matrix for the patch
#' # enriched for CD8
#' j <-
#'   sample(
#'     x = seq_len(length(patch_props)),
#'     size = length(inds_replace),
#'     prob = patch_props,
#'     replace = TRUE
#'   )
#' i <- seq_len(length(inds_replace))
#' x <- rep(1, length(i))
#' dims <- c(length(inds_replace), length(patch_props))
#' dimnames <- list(rownames(counts_with_patch)[inds_replace], names(patch_props))
#'
#' counts <- Matrix::sparseMatrix(
#'   i = i,
#'   j = j,
#'   x = x,
#'   dims = dims,
#'   dimnames = dimnames
#' )
#'
#' counts_with_patch[inds_replace, ] <- counts
#' cg@counts <- counts_with_patch
#'
#'
#' # Run patch detection
#' cg <- patch_detection(
#'   cg,
#'   patch_markers = "CD8"
#' )
#'
#' # Visualize patch
#' xyz <- cg@layout$wpmds_3d %>%
#'   mutate(patch = cg@cellgraph %>% pull(patch))
#'
#' plotly::plot_ly(
#'   xyz,
#'   x = ~x, y = ~y, z = ~z,
#'   color = ~ factor(patch),
#'   type = "scatter3d",
#'   mode = "markers",
#'   colors = c("lightgrey", "red"),
#'   marker = list(
#'     size = 2
#'   )
#' )
#'
#' # Check protein composition
#' gg <- tibble(patch = cg@cellgraph %>% pull(patch) %>% as.factor()) %>%
#'   bind_cols(as.matrix(cg@counts)) %>%
#'   group_by(patch) %>%
#'   summarize(across(where(is.numeric), ~ sum(.x))) %>%
#'   tidyr::pivot_longer(where(is.numeric)) %>%
#'   group_by(patch) %>%
#'   mutate(value = value / sum(value))
#' lvls <- gg %>% arrange(desc(patch), value) %>% pull(name) %>% unique()
#' ggplot(gg %>% mutate(name = factor(name, lvls)), aes(name, patch, fill = value)) +
#'   geom_tile() +
#'   scale_fill_gradientn(
#'     colours = c("lightgrey", "mistyrose", "red", "darkred"),
#'     label = scales::percent
#'   ) +
#'   theme_bw() +
#'   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#'   labs(fill = "Fraction of\ncounts")
#'
#' @return A \code{CellGraph} object with the patches detected as an additional
#' node column in the `tbl_graph` object.
#'
#' @export
#'
patch_detection <- function(
  cg,
  patch_markers,
  receiver_markers = NULL,
  k = 2L,
  leiden_resolution = 0.005,
  patch_nodes_threshold = 100,
  prune_patch_edge = TRUE,
  leiden_refinement = TRUE,
  method = c("expand_contract", "local_G"),
  contraction = 0.5,
  pval_threshold = 0.01,
  seed = 123,
  verbose = TRUE
) {

  # Validate input parameters
  method <- match.arg(method, choices = c("expand_contract", "local_G"))
  assert_class(cg, "CellGraph")
  assert_vector(patch_markers, type = "character", n = 1)
  assert_vector(receiver_markers, type = "character", n = 1, allow_null = TRUE)
  assert_x_in_y(patch_markers, colnames(cg@counts))
  assert_x_in_y(receiver_markers, colnames(cg@counts), allow_null = TRUE)
  assert_single_value(k, type = "integer")
  assert_within_limits(k, limits = c(1L, 4L))
  assert_single_value(leiden_resolution, type = "numeric")
  assert_within_limits(leiden_resolution, limits = c(0, 1))
  assert_single_value(prune_patch_edge, type = "bool")
  assert_single_value(leiden_refinement, type = "bool")
  assert_single_value(contraction, type = "numeric")
  assert_within_limits(contraction, limits = c(0, 1))
  assert_single_value(pval_threshold, type = "numeric")
  assert_within_limits(pval_threshold, limits = c(0, 1))
  assert_single_value(seed, type = "integer")
  assert_single_value(verbose, type = "bool")

  set.seed(seed)

  if (verbose) cli::cli_alert_info(
    "Extracting connected patches using method {.str {method}}..."
  )

  # Select patch detection method
  patch_detection_fkn <- switch(method,
                                expand_contract = .patch_detection_expand_contract,
                                local_G = .patch_detection_local_G)
  s <- switch(method,
              expand_contract = contraction,
              local_G = pval_threshold)

  # Run patch detection
  g_patch <- patch_detection_fkn(cg, patch_markers, receiver_markers, k, prune_patch_edge, s)

  if (verbose) cli::cli_alert(
    c("   {.val {length(g_patch)}} out of {.val {nrow(cg@counts)}}",
      " nodes are labelled as patch nodes...")
  )

  # Split graph into connected components and remove small components with fewer
  # than patch_nodes_threshold nodes
  g_patch <- g_patch %>%
    mutate(comp = igraph::components(.)$membership) %>%
    group_by(comp) %>%
    mutate(n = n())
  # Mark small communities as potential patches
  small_patches <- g_patch %N>%
    as_tibble() %>%
    filter(between(n, 2, patch_nodes_threshold)) %>%
    group_by(comp) %>%
    mutate(potential_patch = cur_group_id()) %>%
    ungroup() %>%
    select(-n, -comp)
  g_exp_patch_list <- g_patch %>%
    filter(n >= patch_nodes_threshold) %>%
    to_components()

  if (length(g_exp_patch_list) == 0) {
    cli::cli_warn(c("i" = "No connected patches found..."))
    cg@cellgraph <- cg@cellgraph %N>%
      mutate(patch = 0L, potential_patch = NA_integer_)
    return(cg)
  }

  if (verbose) cli::cli_alert(
    c("   Found {.val {length(g_exp_patch_list)}}",
      " connected patches after filtering...")
  )

  # Run iterative leiden to split up weakly connected patches
  if (leiden_refinement) {
    if (verbose) cli::cli_alert(
      c("   Splitting up weakly connected patches using",
        " Leiden with resolution={.val {leiden_resolution}}...")
    )
    g_patch_nodes <- lapply(seq_along(g_exp_patch_list), function(i) {
      g_exp_patch_list[[i]] %>%
        mutate(community = igraph::cluster_leiden(
          .,
          resolution = leiden_resolution,
          objective_function = "modularity"
        )$membership) %>%
        as_tibble() %>%
        mutate(patch = paste0(community, "_", i))
    }) %>%
      bind_rows() %>%
      select(-contains("community")) %>%
      group_by(patch) %>%
      mutate(n = n()) %>%
      filter(n >= patch_nodes_threshold)
    if (verbose) cli::cli_alert_info(
      c("   Found {.val {length(unique(g_patch_nodes$patch))}}",
        " patches after splitting...")
    )
  } else {
    g_patch_nodes <- lapply(seq_along(g_exp_patch_list), function(i) {
      g_exp_patch_list[[i]] %N>% as_tibble() %>%
        mutate(patch = paste0("0_", i))
    }) %>%
      bind_rows()
  }

  g_patch_nodes <- g_patch_nodes %>%
    group_by(patch) %>%
    count(name = "size") %>%
    arrange(-size) %>%
    ungroup() %>%
    mutate(group = as.integer(row_number())) %>%
    right_join(g_patch_nodes, by = "patch") %>%
    select(-size) %>%
    ungroup() %>%
    select(name, patch = group)

  # Add group labels to original graph
  g <- cg@cellgraph %N>%
    select(-matches("patch")) %>%
    left_join(g_patch_nodes, by = c("name" = "name")) %>%
    mutate(
      patch = case_when(
        is.na(patch) ~ 0L,
        TRUE ~ patch
      )
    )

  # Add additional small patch node labels
  if (nrow(small_patches) > 0) {
    g <- g %>%
      left_join(small_patches, by = "name")
  }

  if (verbose) cli_alert_success("Finished!")

  cg@cellgraph <- g
  return(cg)
}


#' Detect patches in a CellGraph using the expand-contract method
#'
#' @param cg A `CellGraph` object containing the cell graph.
#' @param patch_markers A character vector of patch protein markers.
#' @param receiver_markers A character vector of host protein markers.
#' @param k The neighborhood size to consider.
#' @param remove_1_core_patch_nodes Logical, whether to remove nodes with core number 1
#' that have multiple connections with the host graph.
#' @param s Scaling factor for the expand-contract method.
#'
#'
#' @noRd
.patch_detection_expand_contract <- function(
  cg,
  patch_markers,
  receiver_markers,
  k,
  remove_1_core_patch_nodes,
  s = 0.5
) {
  # Get adjacency matrix
  A <- cg@cellgraph %>% igraph::as_adjacency_matrix()

  # Expand adjacency matrix
  A_exp <- expand_adjacency_matrix(A, k = k)

  # Define initial patch nodes
  patch_nodes_i <- (cg@counts[, patch_markers, drop = FALSE] %>% Matrix::rowSums()) == 1

  # Define host nodes to block later if desired
  if (!is.null(receiver_markers)) {
    host_nodes <- (cg@counts[, receiver_markers, drop = FALSE] %>% Matrix::rowSums()) == 1
  } else {
    host_nodes <- rep(FALSE, nrow(cg@counts))
  }

  # Update patch_nodes to include outside nodes connecting at least
  # 2 patch nodes within k steps (using expanded adjacency matrix)
  k_step_nodes <- (A_exp[patch_nodes_i, ] %>% Matrix::colSums()) > 1L
  patch_nodes <- patch_nodes_i | k_step_nodes

  # Update patch_nodes to include outside nodes connecting at least
  # 2 patch nodes within 1 steps (using the adjacency matrix)
  k_step_nodes <- Matrix::colSums(A[patch_nodes, ]) > 1L
  patch_nodes <- patch_nodes | k_step_nodes

  # Calculate in/out degree ratio for patch_nodes.
  # Contract: only keep nodes with a ratio >= 1 or the top ranked nodes.
  in_degree <- Matrix::rowSums(A[patch_nodes, patch_nodes])
  out_degree <- Matrix::rowSums(A[patch_nodes, !patch_nodes])
  ratio <- (in_degree + 1) / (out_degree + 1)
  ord_keep <- 1:round(sum(patch_nodes) / (s + 1))
  patch_nodes_keep <- (order(ratio, decreasing = TRUE) %in% ord_keep) | (ratio >= 1)

  if (!is.null(receiver_markers)) {
    host_nodes <- host_nodes &
      (((A[, !patch_nodes] %>% Matrix::rowSums()) > 0) |
      ((A[, host_nodes] %>% Matrix::rowSums()) > 0))
  }

  # Define final patch nodes and subset adjacency matrix
  patch_nodes_final <- (rownames(A) %in% names(patch_nodes_keep)[patch_nodes_keep]) & !host_nodes
  A_patch <- A[patch_nodes_final, patch_nodes_final]

  # Construct initial patch graph
  g_patch <- igraph::graph_from_adjacency_matrix(A_patch, mode = "undirected") %>%
    as_tbl_graph(directed = FALSE)

  # Remove nodes with out degree > 1 and coreness == 1 which are located
  # near the edge(s) of the patch(es). These are more likely to belong to
  # the host graph.
  if (remove_1_core_patch_nodes) {
    g_patch <- .core_1_filter(g_patch, A, patch_nodes_i, patch_nodes_final)
  }

  # Remove orphaned nodes, i.e. nodes that are not connected to any other patch nodes
  A_patch <- igraph::as_adjacency_matrix(g_patch)
  orphaned_patch_nodes <- Matrix::rowSums(A_patch) == 0
  nodes_keep <- names(orphaned_patch_nodes[!orphaned_patch_nodes])
  g_patch <- g_patch %>%
    filter(name %in% nodes_keep)

  return(g_patch)
}


#' Detect patches in a CellGraph using the local G method
#'
#' @param cg A `CellGraph` object containing the cell graph.
#' @param patch_markers A character vector of patch protein markers.
#' @param receiver_markers A character vector of host protein markers.
#' @param k The neighborhood size to consider.
#' @param remove_1_core_patch_nodes Logical, whether to remove nodes with core number 1
#' that have multiple connections with the host graph.
#' @param s Scaling factor for the local G method.
#'
#'
#' @noRd
.patch_detection_local_G <- function(
  cg,
  patch_markers,
  receiver_markers,
  k,
  remove_1_core_patch_nodes,
  s = 0.01
) {

  # Get adjacency matrix
  A <- cg@cellgraph %>% igraph::as_adjacency_matrix()

  # Define initial patch nodes
  patch_nodes_i <- (cg@counts[, patch_markers, drop = FALSE] %>% Matrix::rowSums()) == 1

  # Block host nodes if desired
  if (!is.null(receiver_markers)) {
    host_nodes <- (cg@counts[, receiver_markers, drop = FALSE] %>% Matrix::rowSums()) == 1
  } else {
    host_nodes <- rep(FALSE, nrow(cg@counts))
  }

  patch_counts <- cg@counts[, patch_markers, drop = FALSE] %>% Matrix::rowSums()
  counts <- matrix(patch_counts, ncol = 1, dimnames = list(rownames(cg@counts), "patch_counts")) %>% as("dgCMatrix")

  gi_mat <-
    local_G(
      g = cg@cellgraph,
      counts = counts,
      k = k,
      use_weights = TRUE,
      type = "gstari",
      return_p_vals = TRUE,
      alternative = "greater"
    )
  patch_nodes_final <- gi_mat$gi_p_mat[, 1] < s
  patch_nodes_final <- patch_nodes_final & !host_nodes

  A_patch <- A[patch_nodes_final, patch_nodes_final]

  # Construct patch graph
  g_patch <- igraph::graph_from_adjacency_matrix(A_patch, mode = "undirected") %>%
    as_tbl_graph(directed = FALSE)

  # Remove nodes with out degree > 1 and coreness == 1 which are located
  # near the edge(s) of the patch(es)
  if (remove_1_core_patch_nodes) {
    g_patch <- .core_1_filter(g_patch, A, patch_nodes_i, patch_nodes_final)
  }

  return(g_patch)
}

#' Filter nodes with core number 1 that have multiple connections with the host graph
#'
#' @param g_patch The patch graph.
#' @param A The adjacency matrix of the original cell graph.
#' @param patch_nodes_i Logical vector indicating the initial patch nodes.
#' @param patch_nodes_final Logical vector indicating the final patch nodes.
#'
#' @noRd
.core_1_filter <- function(g_patch, A, patch_nodes_i, patch_nodes_final) {
  nodes_remove <- ((A[patch_nodes_final, !patch_nodes_final] %>% Matrix::rowSums()) > 1) &
    (igraph::coreness(g_patch) == 1) &
    !patch_nodes_i[patch_nodes_final]
  g_patch <- g_patch %>% filter(!nodes_remove)
  return(g_patch)
}



#' Identify population-specific markers for patch detection
#'
#' Patch detection can be applied to identify patches of one cell type
#' on another. The prerequisite for this analysis is that the patch-specific
#' markers are known. Patch detection can be improved further by also leveraging
#' host-specific markers. This function identifies population-specific markers
#' for the patch and host populations, which can then be used for patch detection.
#'
#' @details
#' Patch detection is sensitive to the selection of markers and therefore requires
#' careful selection. The best markers are both high-abundant and specific to a
#' cell population.
#'
#' The method requires a \code{Seurat} object with a metadata column containing
#' the population information, e.g. a column with cell type labels. We then need
#' to specify the host and target populations, where the host population represent
#' the cells on which the patches are expected to be found and the target population
#' represents the cell type from which the patches originated.
#'
#' As a practical example, let's say that our data represent co-cultured T and B cells,
#' and we anticipate that the T cells have patches of B cells on them. Now we face a
#' challenge because the T cell population contains a lot of B cell markers, making it
#' harder to determine what markers are T-cell specific. In other words, the abundance
#' data is mixed.
#'
#' This method attempts to unmix the abundance data using matrix factorization, estimating
#' the composition of the pure host and target populations. The unmixed abundance profiles
#' are then used to label population-specific markers based on the difference in abundance
#' and minimum frequency. The results are summarized in a table, and an optional plot is drawn
#' to help interpret the results.
#'
#' @param object A \code{Seurat} object
#' @param group_by A string specifying the metadata column to group by
#' @param host_population A string specifying the host population name
#' present in the \code{group_by} column
#' @param target_population A string specifying the target population name
#' present in the \code{group_by} column
#' @param abundance_difference A numeric value specifying how many times
#' higher the abundance of a marker should be in one population relative
#' to the other. This is only used to label the markers in the output.
#' @param min_freq A numeric value specifying the minimum frequency of
#' a protein to be labeled in the output table.
#' @param show_plot Logical, whether to show a plot summarizing the results.
#' @param seed An integer seed for reproducibility
#'
#' @return A \code{tbl_df} with the following columns:
#'   - `marker`: the name of the protein.
#'   - `host_unmixed_freq`: the estimated proportion in the host population after unmixing.
#'   - `target_unmixed_freq`: the estimated proportion in the target population after unmixing.
#'   - `host_freq`: the proportion in the host population.
#'   - `target_freq`: the proportion in the target population.
#'   - `label`: a label indicating whether the protein is a marker for the host
#'              or target population. NA values indicate that the protein is unspecific.
#'
#' @export
#'
identify_markers_for_patch_analysis <- function(
  object,
  group_by,
  host_population,
  target_population,
  abundance_difference = 10,
  min_freq = 1e-2,
  show_plot = TRUE,
  seed = 123
) {

  expect_RcppML()
  expect_ggrepel()
  set.seed(seed)

  assert_class(object, "Seurat")
  assert_single_value(group_by, "string")
  assert_single_value(host_population, "string")
  assert_single_value(target_population, "string")
  assert_col_in_data(group_by, object[[]])
  assert_col_class(group_by, object[[]], c("character", "factor"))
  group_vec <- object[[]] %>% pull(all_of(group_by))
  assert_x_in_y(host_population, group_vec)
  assert_x_in_y(target_population, group_vec)
  assert_single_value(abundance_difference, "numeric")
  assert_within_limits(abundance_difference, c(1, 100))
  assert_single_value(min_freq, "numeric")
  assert_within_limits(min_freq, c(0, 1))
  assert_single_value(show_plot, "bool")
  assert_single_value(seed, "integer")

  cells <- colnames(object)[group_vec %in% c(host_population, target_population)]
  counts <- object %>%
    subset(cells = cells) %>%
    JoinLayers() %>%
    LayerData(layer = "counts") %>%
    as.matrix()

  # Get protein proportions for target and host populations
  group_vec <- group_vec[group_vec %in% c(host_population, target_population)]
  target_props <- counts[, group_vec %in% target_population] %>%
    Matrix::rowSums() %>%
    prop.table()
  host_props <- counts[, group_vec %in% host_population] %>%
    Matrix::rowSums() %>%
    prop.table()

  # Run NMF
  nm <- RcppML::nmf(counts, k = 2, verbose = FALSE)
  w <- nm$w %>% prop.table(margin = 2)
  cor_host_w1 <- cor(w[, 1], host_props)
  cor_target_w1 <- cor(w[, 1], target_props)
  if (max(c(cor_host_w1, cor_target_w1)) < 0.9) {
    cli::cli_abort(
      c(
        "x" = "Failed to identify population markers."
      )
    )
  }
  if (cor_target_w1 > cor_host_w1) {
    target_unmixed_props <- w[, 1]
    host_unmixed_props <- w[, 2]
  } else {
    target_unmixed_props <- w[, 2]
    host_unmixed_props <- w[, 1]
  }

  # Create tibble
  df <- tibble(
    marker = rownames(object),
    tup = target_unmixed_props,
    hup = host_unmixed_props,
    tp = target_props,
    hp = host_props
  ) %>%
    mutate(
      label = case_when(
        (tup > abundance_difference * hup) & (hup < min_freq) ~ glue::glue("marker for {target_population}"),
        (hup > abundance_difference * tup) & (tup < min_freq)  ~ glue::glue("marker for {host_population}"),
        TRUE ~ "Unspecific"
      )
    ) %>%
    rename(
      host_freq = hp,
      target_freq = tp,
      host_unmixed_freq = hup,
      target_unmixed_freq = tup
    )

  # Show plot
  if (show_plot) {
    visualize_patch_analysis_markers(
      df,
      host_population,
      target_population
    ) %>%
      print()
  }

  return(df)
}


#' Visualize patch analysis markers
#'
#' @param df A data frame containing labelled markers
#'
#' @noRd
visualize_patch_analysis_markers <- function(
  df,
  host_population,
  target_population
) {
  dfm <- bind_rows(
    df %>%
      select(marker, host_freq, target_freq, label) %>%
      mutate(type = "True"),
    df %>%
      select(marker, host_unmixed_freq, target_unmixed_freq, label) %>%
      rename(host_freq = host_unmixed_freq, target_freq = target_unmixed_freq) %>%
      mutate(type = "Estimated")
  ) %>%
    mutate(type = factor(type, c("True", "Estimated")))
  cols <- set_names(
    c("#4A73C0", "#DB5D61", "lightgrey"),
    c(glue::glue("marker for {host_population}"),
      glue::glue("marker for {target_population}"),
      "Unspecific")

  )
  p1 <- ggplot(dfm, aes(host_freq, target_freq, color = label, size = pmax(host_freq, target_freq))) +
    geom_point() +
    theme_bw() +
    ggrepel::geom_text_repel(aes(label = marker), max.overlaps = 50) +
    scale_x_continuous(expand = c(0, 0.02), labels = scales::percent) +
    scale_y_continuous(expand = c(0, 0.02), labels = scales::percent) +
    scale_size(range = c(1, 4)) +
    labs(x = glue::glue("Proportion in {host_population}"),
         y = glue::glue("Proportion in {target_population}"),
         color = "") +
    facet_grid(~type) +
    guides(size = "none") +
    scale_color_manual(values = cols)
  est_props <- df %>%
    na.omit() %>%
    group_by(label) %>%
    summarize(hp_tot = sum(host_freq),
              tp_tot = sum(target_freq),
              hup_tot = sum(host_unmixed_freq),
              tup_tot = sum(target_unmixed_freq)) %>%
    pivot_longer(all_of(c("hp_tot", "tp_tot", "hup_tot", "tup_tot"))) %>%
    mutate(type = if_else(
      stringr::str_detect(name, "hup|tup"), "Estimated", "True"
    )) %>%
    mutate(name = case_when(
      name %in% c("hp_tot", "hup_tot") ~ host_population,
      name %in% c("tp_tot", "tup_tot") ~ target_population
    )) %>%
    mutate(type = factor(type, c("True", "Estimated")))
  p2 <- ggplot(est_props, aes(name, value, fill = label)) +
    geom_col() +
    theme_bw() +
    scale_y_continuous(expand = c(0, 0.02), labels = scales::percent) +
    labs(x = "", y = "Proportion of total abundance", fill = "") +
    facet_grid(~type) +
    coord_flip() +
    scale_fill_manual(values = cols)
  wrap_plots(p1, p2, ncol = 1, heights = c(2, 1))
}


#' Expand an adjacency matrix to include higher-order neighborhoods
#'
#' The node neighborhoods are expanded by multiplying the adjacency matrix \code{k}
#' times. To make sure that all neighbors up to \code{k} steps away are included,
#' self-loops are added by setting the diagonal of \code{A} to 1. These self-loops
#' are removed from the resulting matrix after expansion. Note that the self-loops
#' affect how the transition probabilities are calculated. When \code{use_weights = FALSE},
#' the edge values of the result from the matrix multiplication corresponds to the
#' number of k-step random walks between any pair of nodes. For the expanded adjacency
#' matrix, we replace these edge values with 1 so that the adjacency matrix corresponds
#' to an unweighted graph with expanded neighborhoods.
#'
#' @param A An adjacency matrix. Can be a sparse matrix (\code{dgCMatrix})
#' @param k An integer specifying the neighborhood size. Each node will be
#' connected to nodes up to \code{k} steps away.
#' @param use_weights If \code{TRUE}, the transpose of the stochastic matrix
#' will be expanded. The resulting edge weights will correspond to the incoming
#' transition probabilities for a \code{k}-step random walk.
#'
#' @returns An expanded adjacency matrix
#'
#' @export
#'
expand_adjacency_matrix <- function(
  A,
  k = 1L,
  use_weights = FALSE
) {
  assert_class(A, c("matrix", "dgCMatrix"))
  assert_single_value(k, "integer")
  assert_single_value(use_weights, "bool")

  if (!use_weights && k == 1L) {
    return(A)
  }

  if (use_weights && k == 1L) {
    A <- (A / Matrix::rowSums(A)) %>% Matrix::t()
  }

  # Expand neighborhood when k > 1
  if (k > 1L) {
    Matrix::diag(A) <- 1
    if (use_weights) {
      # Use incoming transition probabilities as weights
      A <- (A / Matrix::rowSums(A)) %>% Matrix::t()
    }
    # Expand neighborhood
    A <- Reduce("%*%", rep(list(A), k))
    if (!use_weights) {
      # Set all edge weights to 1 if not using weights
      A@x <- rep(1, length(A@x))
    }
    Matrix::diag(A) <- 0
  }
  return(A)
}


#' Calculate cosine distances between two sets of coordinates
#'
#' @param A,B Matrices of coordinates. Each row corresponds to a point in space.
#'
#' @export
#'
cos_dist <- function(A, B) {
  assert_class(A, "matrix")
  assert_class(B, "matrix")
  if (!ncol(A) == ncol(B)) {
    cli::cli_abort(
      c("i" = "A and B must have the same number of columns",
        x = "{.var A} has {ncol(A)} columns, but {.var B} has {ncol(B)} columns.")
    )
  }

  # Compute dot product
  dot_product <- tcrossprod(A, B)

  # Compute magnitudes
  magnitude_A <- sqrt(rowSums(A^2))
  magnitude_B <- sqrt(rowSums(B^2))

  # Compute cosine distance
  cosine_dist <- (dot_product / (magnitude_A %*% t(magnitude_B)))

  return(cosine_dist)
}
