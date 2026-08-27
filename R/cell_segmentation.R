#' cell:cell conjugate segmentation
#'
#' This function performs segmentation of a cell:cell conjugate graph into its respective
#' cell types and the interface between them (optional).
#'
#' Requirements:
#' - A `CellGraph` object containing the graph representation of the cell:cell conjugate. Only
#' two cell types should be present in the graph.
#' - A matrix of NMF weights (`w`) for the two cell types, which can be derived using the
#' `cc_protein_weights` function.
#'
#' The segmentation is performed in several steps:
#'
#' 1. Classify nodes into "cell1" and "cell2" using Non-Negative Least Squares (NNLS). This
#' method is described in detail below. After classification, the Largest Connected Component (LCC) is
#' kept for each cell type to remove small disconnected subgraphs.
#' 2. (Optional) Identify the interface nodes between the two cell types. These are defined as nodes
#' that are directly connected to both cell types. Optionally, the interface can be expanded
#' by including neighboring nodes within `k_interface_expansion` hops. By setting `detect_interface`
#' to `FALSE`, the interface detection step is skipped.
#'
#' Nodes that cannot be classified with confidence are labeled as "other".
#'
#' @section NNLS method for classifying membrane nodes:
#' To segment cell:cell conjugates, we first need to fit an NMF model to the abundance matrix (or local node
#' neighborhood abundance profiles) of the two cell populations using the `cc_protein_weights` function.
#' This gives us the protein weights for each cell type (`w`) that we can use for classification. As a rule of thumb,
#' a good NMF model should have high weights for proteins that are mutually exclusive between the two cell types.
#' For instance, if we are segmenting a B cell:T cell conjugate, we would expect the NMF model to give high weights
#' to CD10, CD20 and CD21 for the B cell identity and CD3e, CD2 and CD5 for the T cell identity.
#' For the membrane segmentation step, we create a neighborhood profile for each node by summing the abundance
#' of its neighbors in the graph. We then project the NMF weights (`w`) onto these neighborhood profiles using
#' NNLS to get a score for each node and cell type.
#'
#' @param cg A `CellGraph` object
#' @param w A matrix of NMF weights with rows corresponding to proteins and columns
#' corresponding to the two interacting cell types.
#' @param k The number of hops to use when computing local node abundance profiles used to
#' compute the NNLS projection scores.
#' @param detect_interface Whether to identify the interface between the two cell types.
#' @param k_interface_expansion The number of hops to expand the interface between the two cell types.
#' The default (0) classifies nodes that are directly connected to both cell types as interface.
#' Values greater than 0 can be set to expand the interface further by including nodes that are within
#' k hops of the initial interface nodes.
#' @param keep_largest_comp Whether to keep only the largest connected component in each
#' segmented cell-type graph. If `FALSE`, all connected components with at least
#' `min_comp_size` nodes are kept.
#' @param min_comp_size The minimum size of connected components to keep when
#' `keep_largest_comp` is `FALSE`.
#' @param spatial_smoothing_iter The number of iterations to perform for spatial smoothing of
#' node neighborhoods scores. This helps smoothen out local noise in the graph and improve classification.
#' @param verbose Print updates during the segmentation process.
#'
#' @examples
#' library(dplyr)
#' se <- ReadPNA_Seurat(minimal_pna_pxl_file()) %>%
#'   LoadCellGraphs(cells = colnames(.)[2], verbose = FALSE)
#' cg <- CellGraphs(se)[[2]]
#'
#' # Get cell type weights
#' se$cell_type <- c("Mono", "pDC", "CD4T", "CD4T", "CD4T")
#' w <- cc_protein_weights(
#'   se,
#'   group_by = "cell_type",
#'   population_1 = "Mono", population_2 = "CD4T",
#'   show_plot = FALSE
#' )
#'
#' cg <- segment_cell(
#'   cg,
#'   w = w
#' )
#'
#' @return A `CellGraph` object containing the compartment labels from the segmentation as a node attribute.
#'
#' @export
#'
segment_cell <- function(
  cg,
  w,
  k = 2L,
  detect_interface = TRUE,
  k_interface_expansion = 0L,
  keep_largest_comp = TRUE,
  min_comp_size = 10,
  spatial_smoothing_iter = 5L,
  verbose = TRUE
) {
  # Validate input parameters
  .validate_segment_cell_input_params(
    cg,
    w,
    k,
    detect_interface,
    k_interface_expansion,
    keep_largest_comp,
    min_comp_size,
    spatial_smoothing_iter,
    verbose
  )

  if (!detect_interface) {
    k_interface_expansion <- -1
  }

  # Use the column names of w as cell names.
  cell_names <- colnames(w)

  # Fetch the adjacency matrix of the graph
  A <- igraph::as_adjacency_matrix(cg@cellgraph)

  # Expand the adjacency matrix to create a k-hop reachability matrix.
  A_k <- expand_adjacency_matrix(A, k = k)
  Matrix::diag(A_k) <- 1

  # Use all nodes for scoring; component filtering is applied after classification.
  counts_segment <- cg@counts
  A_k <- A_k[rownames(counts_segment), rownames(counts_segment)]

  # Classify nodes based on NNLS projection scores.
  node_classification <- .classify_nodes_nnls(
    counts_membrane = as.matrix(counts_segment),
    A_k = A_k,
    A = A[rownames(A_k), colnames(A_k)],
    cell_names = cell_names,
    w = w,
    spatial_smoothing_iter = spatial_smoothing_iter,
    verbose = verbose
  )

  # Create cell1 and cell2 graphs by subsetting on classification and filtering components.
  # Index into node_classification by name to guard against any ordering divergence
  # between cg@counts rows and cg@cellgraph nodes.
  graph_node_names <- cg@cellgraph %N>% pull(name)
  g_c1 <- .subset_graph(
    cg@cellgraph,
    bool_filter = node_classification[graph_node_names] == 1,
    keep_largest_comp = keep_largest_comp,
    min_comp_size = min_comp_size
  )
  g_c2 <- .subset_graph(
    cg@cellgraph,
    bool_filter = node_classification[graph_node_names] == 2,
    keep_largest_comp = keep_largest_comp,
    min_comp_size = min_comp_size
  )

  if (verbose) {
    cli::cli_alert_info(
      "Defining interface nodes between {.str {cell_names}}"
    )
  }

  if (k_interface_expansion == -1) {
    if (verbose) {
      cli::cli_alert_info(
        c(
          "Skipping interface detection."
        )
      )
    }
    interface_nodes <- NULL
  } else {
    interface_nodes <- .fetch_interface_nodes(
      A,
      c1_nodes = g_c1 %>% pull(name),
      c2_nodes = g_c2 %>% pull(name),
      k_interface_expansion = k_interface_expansion
    )
  }

  # Map nodes to compartments
  # NOTE: currently, the interface breaks the cell1 and cell2 graphs
  if (verbose) {
    cli::cli_alert_info(
      "Mapping nodes to compartments"
    )
  }

  # Final mapping. The priority of the mapping is as follows:
  # 1) Interface nodes are classified as "interface"
  # 2) Nodes in retained cell-specific components are classified as "cell1" or "cell2"
  # 3) Nodes that don't fit any of the above categories are classified as "other"
  # Note that the interface can break the connectivity of the cell1 and cell2 graphs,
  # so there's no guarantee that these compartments will be fully connected.
  c1_nodes_keep <- g_c1 %>% pull(name)
  c2_nodes_keep <- g_c2 %>% pull(name)

  node_compartment_map <- tibble(
    node = rownames(cg@counts)
  ) %>%
    mutate(group = if_else(
      node %in% interface_nodes, "interface", "other"
    )) %>%
    mutate(compartment = case_when(
      group == "interface" ~ "interface",
      node %in% c1_nodes_keep ~ cell_names[1],
      node %in% c2_nodes_keep ~ cell_names[2],
      TRUE ~ "other"
    ))

  # Add compartment annotation to graph
  cg@cellgraph <- cg@cellgraph %N>%
    select(-any_of("compartment")) %>%
    left_join(node_compartment_map %>% select(node, compartment), by = c("name" = "node"))

  return(cg)
}


#' Validate input parameters for segment_cell function
#'
#' @param cg A `CellGraph` object
#' @param w A matrix of NMF weights with rows corresponding to proteins and columns
#' corresponding to the two interacting cell types.
#' @param k The number of hops to use when computing cell type identity scores.
#' @param detect_interface Whether to identify the interface between the two cell types.
#' @param k_interface_expansion The number of hops to expand the interface between the two cell types.
#' @param keep_largest_comp Whether to keep only the largest connected component in each segmented cell-type graph.
#' @param min_comp_size The minimum size of connected components to keep when `keep_largest_comp` is `FALSE`.
#' @param spatial_smoothing_iter The number of iterations to perform for spatial smoothing of node neighborhoods scores.
#' @param verbose Print updates during the segmentation process.
#' @param call The calling environment for error messages. Defaults to the caller's environment.
#'
#' @return Nothing. Only used for its side effects
#'
#' @noRd
.validate_segment_cell_input_params <- function(
  cg,
  w,
  k,
  detect_interface,
  k_interface_expansion,
  keep_largest_comp,
  min_comp_size,
  spatial_smoothing_iter,
  verbose,
  call = rlang::caller_env()
) {
  assert_class(cg, "CellGraph", call = call)
  if (is.null(w)) {
    cli::cli_abort(
      c(
        "x" = "{.var w} must be provided for cell:cell segmentation."
      ),
      call = call
    )
  }
  assert_class(w, "matrix", call = call)
  assert_vector(colnames(w), type = "character", n = 2, call = call)
  assert_vector(rownames(w), type = "character", n = 1, call = call)
  if (ncol(w) != 2L) {
    cli::cli_abort(
      c(
        "x" = "{.var w} must have exactly two columns corresponding to the two interacting cell types"
      ),
      call = call
    )
  }
  assert_single_value(k, type = "integer", call = call)
  assert_within_limits(k, limits = c(1, 6), call = call)
  assert_single_value(detect_interface, type = "bool", call = call)
  assert_single_value(k_interface_expansion, type = "integer", call = call)
  assert_within_limits(k_interface_expansion, limits = c(0, 4), call = call)
  assert_single_value(keep_largest_comp, type = "bool", call = call)
  assert_single_value(min_comp_size, type = "integer", call = call)
  assert_within_limits(min_comp_size, limits = c(1, Inf), call = call)
  assert_single_value(spatial_smoothing_iter, type = "integer", call = call)
  assert_within_limits(spatial_smoothing_iter, limits = c(0, 50), call = call)
  assert_single_value(verbose, type = "bool", call = call)
}

#' Derive protein weights for two cell populations
#'
#' This function uses non-negative matrix factorization (NMF) to derive protein
#' weights for two cell populations based on their protein abundance. It identifies
#' proteins that are differentially expressed between the two populations and assigns
#' weights that can be used for cell type inference.
#'
#' By default, a plot is generated showing the top protein weights for each population,
#' which can help evaluate the specificity of the identified markers and cross-reference
#' these with literature.
#'
#' @section "k_neighborhood" mode details:
#' In "k_neighborhood" mode, the NMF model is fit on count profiles from local graph
#' neighborhoods instead of whole-cell abundance profiles. Because `segment_cell()`
#' classifies nodes from neighborhood-level information, this mode can improve marker
#' specificity and downstream segmentation quality.
#'
#' This mode is more computationally expensive and typically requires tuning of
#' neighborhood-related parameters (for example `k`, sampling depth, and minimum
#' neighborhood size).
#'
#' Both "cell_abundance" and "k_neighborhood" modes assume that cell-type abundance
#' differences are the dominant source of variation. If other effects are stronger
#' (for example spatial structure), NMF may capture those effects instead of cell-type
#' signal. This risk is often higher in "k_neighborhood" mode because local profiles
#' are more sensitive to spatial patterns.
#'
#' To reduce this risk, increase `k` (to smooth local variation) and/or provide
#' `masked_markers` to exclude known confounding markers from model fitting.
#'
#' When using "k_neighborhood" mode, neighborhoods are split into a training set (80%)
#' and a test set (20%). The NMF model is fit on the training set, and the learned
#' weights are then applied to the test set to assess generalization. Low test-set
#' classification performance suggests the model may not be capturing stable,
#' cell type-specific signal.
#'
#' @section masking markers:
#' If strong spatial signatures interfere with cell-type separation, supply
#' `masked_markers` to exclude those proteins during model fitting. A common example
#' is platelet contamination, where markers such as CD41, CD36, CD62P, and CD9 are
#' often masked.
#'
#' @param object A `Seurat` object with cell type labels in the metadata.
#' @param group_by The name of the metadata column containing the cell type
#' population labels.
#' @param population_1,population_2 The names of the two cell populations.
#' @param masked_markers Optional character vector of proteins to exclude from the model.
#' @param max_freq Maximum frequency of a protein in either population to be
#' included in the model. High abundant proteins such as HLA-ABC tend to get
#' high protein weights across both cell populations, making them unsuitable
#' for cell type inference. Default is 0.01 (1%).
#' @param min_diff Minimum difference in protein frequency between the two
#' populations. This is calculated as (freq_pop1 - freq_pop2) / (freq_pop1 + freq_pop2).
#' @param mode Method to use for deriving protein weights. "k_neighborhood" uses
#' the abundance of proteins in k-hop neighborhoods, while "cell_abundance" uses
#' the overall abundance of proteins in each cell population. Default is "cell_abundance".
#' @param max_components_per_population Maximum number of components (cells or neighborhoods)
#' to sample from each population for model training when `mode` is "k_neighborhood". Default is 100.
#' @param k The neighborhood size to use when `mode` is "k_neighborhood". This
#' defines the number of hops to consider when creating the neighborhood profiles
#' for each node. Default is 2.
#' @param neighborhoods_per_component The number of neighborhoods to sample per cell
#' type when `mode` is "k_neighborhood".
#' @param min_neighborhood_size Minimum number of nodes in a neighborhood for it
#' to be included in the model when `mode` is "k_neighborhood". Default is 10.
#' @param verbose Whether to print progress messages during the weight derivation process.
#'
#' @param show_plot Whether to show a plot of the top protein weights.
#' Default is TRUE.
#'
#' @examples
#' library(dplyr)
#' se <- ReadPNA_Seurat(minimal_pna_pxl_file())
#' se$cell_type <- c("Mono", "pDC", "CD4T", "CD4T", "CD4T")
#' w <- cc_protein_weights(
#'   se,
#'   group_by = "cell_type",
#'   population_1 = "Mono",
#'   population_2 = "CD4T"
#' )
#'
#' head(w)
#'
#' @return A matrix of protein weights for the two populations.
#'
#' @export
#'
cc_protein_weights <- function(
  object,
  group_by,
  population_1,
  population_2,
  masked_markers = NULL,
  max_freq = 0.01,
  min_diff = 0.5,
  mode = c("cell_abundance", "k_neighborhood"),
  max_components_per_population = 100L,
  neighborhoods_per_component = 500L,
  min_neighborhood_size = 10L,
  k = 2L,
  show_plot = TRUE,
  verbose = TRUE
) {
  expect_RcppML()

  .validate_cc_protein_weights_input_params(
    object = object,
    group_by = group_by,
    population_1 = population_1,
    population_2 = population_2,
    masked_markers = masked_markers,
    max_freq = max_freq,
    min_diff = min_diff,
    max_components_per_population = max_components_per_population,
    neighborhoods_per_component = neighborhoods_per_component,
    min_neighborhood_size = min_neighborhood_size,
    k = k,
    show_plot = show_plot,
    verbose = verbose
  )

  mode <- match.arg(mode, c("cell_abundance", "k_neighborhood"))

  # Define group vector and check that specified populations are present in the data
  group_vec <- object[[]] %>% pull(!!sym(group_by))
  assert_x_in_y(population_1, group_vec)
  assert_x_in_y(population_2, group_vec)

  # List of length 2 with logical vectors indicating which components (cells or neighborhoods)
  # to keep for each population.
  pop_components <- .define_components_to_keep(
    group_vec,
    population_1,
    population_2,
    max_components_per_population,
    verbose
  )

  # Build count matrix from local neighborhoods
  if (mode == "k_neighborhood") {
    results <- .build_count_matrix_for_k_neighborhood(
      object,
      pop_components,
      masked_markers,
      k,
      neighborhoods_per_component,
      min_neighborhood_size,
      max_components_per_population,
      verbose
    )
  }

  # Build count matrix from the cell abundance matrix
  if (mode == "cell_abundance") {
    results <- .build_count_matrix_for_cell_abundance(
      object,
      pop_components
    )
  }

  counts <- results$counts
  pop_components <- results$pop_components

  # Compute marker frequencies for populations
  pop1_props <- counts[, pop_components[[1]], drop = FALSE] %>%
    Matrix::rowSums() %>%
    prop.table()
  pop2_props <- counts[, pop_components[[2]], drop = FALSE] %>%
    Matrix::rowSums() %>%
    prop.table()

  # Filter out unspecific markers to improve model
  diff <- (pop1_props - pop2_props) / (pop1_props + pop2_props)
  diff[!is.finite(diff)] <- 0
  markers_keep <- abs(diff) > min_diff
  markers_keep <- markers_keep & (pmin(pop1_props, pop2_props) < max_freq)
  if (!is.null(masked_markers)) {
    markers_keep <- markers_keep & !(rownames(counts) %in% masked_markers)
  }

  # Check how many markers are kept after filtering
  if (sum(markers_keep) < 5) {
    cli::cli_abort(
      c(
        "x" = "{sum(markers_keep)} markers are kept after filtering with {.var max_freq} and {.var min_diff}. ",
        "i" = "Consider relaxing the filtering criteria and check that the populations are correctly specified.",
        "!" = "Note that segmentation is not an option for similar cell types."
      )
    )
  }
  if (sum(markers_keep) < 10) {
    cli::cli_alert_warning(
      c(
        "i" = "{sum(markers_keep)} markers are kept after filtering with {.var max_freq} and {.var min_diff}. ",
        "!" = "This is a low number of markers and may lead to unstable results. "
      )
    )
  }

  # Filter count matrix
  counts <- as.matrix(counts[markers_keep, ])

  if (mode == "k_neighborhood") {
    # Only keep neighborhoods that have at least min_neighborhood_size nodes after
    # filtering to ensure robust model training.
    neighborhoods_to_keep <- Matrix::colSums(counts) >= min_neighborhood_size
    counts <- counts[, neighborhoods_to_keep, drop = FALSE]
    pop_components <- lapply(pop_components, function(x) x[neighborhoods_to_keep])

    # Ensure there are enough neighborhoods remaining after filtering
    if (ncol(counts) < 1e3) {
      cli::cli_abort(
        paste0(
          "After filtering, too few neighborhoods ({ncol(counts)}",
          ") remain to train and evaluate an NMF model in 'k_neighborhood' mode. ",
          "Consider increasing {.var k} or use mode 'cell_abundance'."
        )
      )
    }

    # Split into training and test sets to evaluate model performance. We use 80%
    # of the data for training and 20% for testing.
    train_indices <- sample(seq_len(ncol(counts)), size = 0.8 * ncol(counts), replace = FALSE)
    test_indices <- setdiff(seq_len(ncol(counts)), train_indices)
    if (verbose) {
      cli::cli_alert_info(
        c("After filtering, {ncol(counts)} neighborhoods remain for model training and testing.")
      )
      cli::cli_alert_info(
        c(
          "Splitting data into training (80%) and test (20%) sets with ",
          "{length(train_indices)} and {length(test_indices)} neighborhoods, respectively."
        )
      )
    }
    counts_train <- counts[, train_indices, drop = FALSE]
    counts_test <- counts[, test_indices, drop = FALSE]
    pop_components_test <- lapply(pop_components, function(x) x[test_indices])
    pop_components <- lapply(pop_components, function(x) x[train_indices])

    # Fit model on training data
    nmf_model <- RcppML::nmf(counts_train, k = 2, L1 = c(0, 0), tol = 1e-10, verbose = FALSE)
    w <- nmf_model$w
  }
  if (mode == "cell_abundance") {
    # Fit model
    nmf_model <- RcppML::nmf(counts, k = 2, L1 = c(0, 0), tol = 1e-10, verbose = FALSE)
    w <- nmf_model$w
  }

  # Ensure correct ordering of weight vectors
  pop1_ind <- nmf_model$h[, pop_components[[1]], drop = FALSE] %>%
    apply(1, sum) %>%
    which.max()
  if (pop1_ind == 2) {
    w <- w[, 2:1]
  }
  colnames(w) <- c(population_1, population_2)
  rownames(w) <- rownames(counts)

  if (show_plot) {
    gg <- tibble(
      score = c(w[, population_1], w[, population_2]),
      marker = rep(rownames(w), 2),
      component = rep(c(population_1, population_2), each = nrow(w))
    )

    top_markers <- gg %>%
      group_by(component) %>%
      slice_max(score, n = 40) %>%
      pull(marker) %>%
      unique()

    p <- gg %>%
      filter(marker %in% top_markers) %>%
      ggplot(aes(reorder(marker, score), score, fill = component)) +
      geom_col() +
      coord_flip() +
      theme_bw() +
      labs(fill = "Population", x = "Protein", y = "Weights", title = "NMF model weights") +
      scale_fill_manual(values = c("#6588CF", "#F19DA7") %>% set_names(c(population_1, population_2)))

    print(p)
  }

  if (mode == "k_neighborhood") {
    if (verbose) {
      cli::cli_alert_info(
        "Evaluating model performance on test set"
      )
    }
    # Evaluate model performance
    x <- .score_cell_types(
      counts = counts_test,
      w = w
    )

    km <- kmeans(x[, population_1], centers = 2)
    th <- .kmeans_midpoint(km$centers)
    classification <- if_else(x[, population_1] >= th, population_1, population_2) %>%
      set_names(rownames(x))

    classification_table <- tibble(
      actual = pop_components_test[[1]] %>% if_else(., population_1, population_2),
      predicted = classification
    ) %>%
      group_by(actual, predicted) %>%
      count() %>%
      group_by(actual) %>%
      mutate(p = n / sum(n)) %>%
      mutate(type = factor(if_else(actual == predicted, "accuracy", "error"), levels = c("accuracy", "error"))) %>%
      select(actual, type, p) %>%
      pivot_wider(names_from = type, values_from = p, values_fill = list(p = 0), names_expand = TRUE) %>%
      ungroup() %>%
      data.frame(row.names = 1, check.names = FALSE) %>%
      as.matrix()

    if (verbose) {
      cli::cli_alert_info(
        c(
          "Model accuracy on test set: ",
          "{classification_table[population_1, 'accuracy'] %>% scales::percent(1)} for {.str {population_1}} and ",
          "{classification_table[population_2, 'accuracy'] %>% scales::percent(1)} for {.str {population_2}}."
        )
      )
    }

    # Store model type and parameters as attributes of w
    attributes(w)$model_type <- "k_neighborhood"
    attributes(w)$k <- k
    attributes(w)$neighborhoods_per_component <- neighborhoods_per_component
    attributes(w)$min_neighborhood_size <- min_neighborhood_size
    attributes(w)$classification_table <- classification_table
  } else {
    attributes(w)$model_type <- "cell_abundance"
  }

  return(w)
}

#' Internal function to define components to keep for model training
#'
#' @return A list with two logical vectors indicating which components (cells) to keep for each population.
#'
#' @noRd
.define_components_to_keep <- function(
  group_vec,
  population_1,
  population_2,
  max_components_per_population,
  verbose
) {
  components_keep_pop_1 <- group_vec == population_1
  components_keep_pop_2 <- group_vec == population_2

  n_pop_1 <- sum(components_keep_pop_1)
  n_pop_2 <- sum(components_keep_pop_2)
  if (verbose) {
    if (n_pop_1 < 20) {
      cli::cli_alert_warning(
        c("Population {.str {population_1}} has less than 20 cells, which may lead to unstable results.")
      )
    }
    if (n_pop_2 < 20) {
      cli::cli_alert_warning(
        c("Population {.str {population_2}} has less than 20 cells, which may lead to unstable results.")
      )
    }
  }
  if (n_pop_1 > max_components_per_population) {
    inds_keep1 <- sample(which(components_keep_pop_1), size = max_components_per_population)
    components_keep_pop_1 <- rep(FALSE, length(components_keep_pop_1))
    components_keep_pop_1[inds_keep1] <- TRUE
  }
  if (n_pop_2 > max_components_per_population) {
    inds_keep2 <- sample(which(components_keep_pop_2), size = max_components_per_population)
    components_keep_pop_2 <- rep(FALSE, length(components_keep_pop_2))
    components_keep_pop_2[inds_keep2] <- TRUE
  }

  return(list(components_keep_pop_1, components_keep_pop_2) %>% set_names(c(population_1, population_2)))
}

#' Internal function to build count matrix for cell abundance mode of `cc_protein_weights`
#'
#' @return A count matrix
#'
#' @noRd
.build_count_matrix_for_cell_abundance <- function(
  object,
  pop_components
) {
  assay <- DefaultAssay(object)
  pxl_assay <- object[[assay]]

  if (inherits(pxl_assay, "PNAAssay5")) {
    if (length(Layers(pxl_assay, layer = "counts") %>% grep("^counts\\.\\d", x = .)) > 1) {
      cli::cli_abort(
        c(
          "Multiple count layers found in {.var object}. Please run {.var JoinLayers} and try again."
        )
      )
    }
  }
  # Fetch count matrix
  pop1_components <- colnames(object)[pop_components[[1]]]
  pop2_components <- colnames(object)[pop_components[[2]]]
  counts <- LayerData(object, layer = "counts")[, c(pop1_components, pop2_components), drop = FALSE]

  pop_components <- list(
    c(rep(TRUE, length(pop1_components)), rep(FALSE, length(pop2_components))),
    c(rep(FALSE, length(pop1_components)), rep(TRUE, length(pop2_components)))
  )

  return(list(
    counts = counts,
    pop_components = pop_components
  ))
}

#' Internal function to build count matrix for k-neighborhood mode of `cc_protein_weights`
#'
#' @return A list with:
#' - `counts` a count matrix
#' - `group_vec` an updated group vector
#'
#' @noRd
.build_count_matrix_for_k_neighborhood <- function(
  object,
  pop_components,
  masked_markers,
  k,
  neighborhoods_per_component,
  min_neighborhood_size,
  max_components_per_population,
  verbose
) {
  if (verbose) {
    cli::cli_alert_info(
      c(
        "Fetching {neighborhoods_per_component} {k}-step neighborhoods from ",
        "{sum(pop_components[[1]])} {.str {names(pop_components)[1]}} cells and ",
        "{sum(pop_components[[2]])} {.str {names(pop_components)[2]}} cells."
      )
    )
  }

  pop1_components <- colnames(object)[pop_components[[1]]]
  pop2_components <- colnames(object)[pop_components[[2]]]

  # Fetch neighborhoods separately per population so that population membership
  # is derived from the actual returned columns, not positional arithmetic.
  # SQL LIMIT means cells with fewer nodes than neighborhoods_per_component return
  # fewer seeds; a joint call followed by a fixed-size positional split would
  # silently misalign the population masks.
  counts_pop1 <- fetch_node_neighborhoods(
    object,
    pop1_components,
    nodes_per_component = neighborhoods_per_component,
    k = k
  )
  counts_pop2 <- fetch_node_neighborhoods(
    object,
    pop2_components,
    nodes_per_component = neighborhoods_per_component,
    k = k
  )

  counts <- SeuratObject::RowMergeSparseMatrices(counts_pop1, list(counts_pop2))

  # Build population masks from actual column names, not expected counts
  pop_components <- list(
    colnames(counts) %in% colnames(counts_pop1),
    colnames(counts) %in% colnames(counts_pop2)
  ) %>% set_names(names(pop_components))

  # Check percent masked
  if (!is.null(masked_markers)) {
    pct_masked <- Matrix::colSums(counts[
      intersect(masked_markers, rownames(counts)), ,
      drop = FALSE
    ]) / Matrix::colSums(counts)
    nbs_keep_clean <- pct_masked < 0.5
    if ((sum(!nbs_keep_clean) / length(nbs_keep_clean)) > 0.1) {
      cli::cli_alert_warning(
        c(
          "More than 10% of neighborhoods have high `masked_markers` content. ",
          "This could lead to unstable results."
        )
      )
    }
  }
  # Filter by minimum allowed neighborhood size
  nbs_keep_large <- Matrix::colSums(counts) >= min_neighborhood_size
  if ((sum(!nbs_keep_large) / length(nbs_keep_large)) >= 0.1) {
    cli::cli_alert_warning(
      c(
        "More than 10% of neighborhoods will be discarded from the model training ",
        "step due to small neighborhood sizes. Consider increasing {.var k} to get larger neighborhoods."
      )
    )
  }

  nbs_keep <- nbs_keep_large
  if (!is.null(masked_markers)) {
    nbs_keep <- nbs_keep & nbs_keep_clean
  }

  counts <- counts[, nbs_keep, drop = FALSE]
  pop_components <- lapply(pop_components, function(x) x[nbs_keep])

  return(list(
    counts = counts,
    pop_components = pop_components
  ))
}


#' Internal function to validate input parameters for `cc_protein_weights`
#'
#' @return Nothing. Only used for validation.
#'
#' @noRd
.validate_cc_protein_weights_input_params <- function(
  object,
  group_by,
  population_1,
  population_2,
  masked_markers,
  max_freq,
  min_diff,
  show_plot,
  k,
  min_neighborhood_size,
  max_components_per_population,
  neighborhoods_per_component,
  verbose,
  call = rlang::caller_env()
) {
  assert_class(object, "Seurat", call = call)
  assert_single_value(group_by, type = "string", call = call)
  assert_col_in_data(group_by, object[[]], call = call)
  assert_single_value(population_1, type = "string", call = call)
  assert_single_value(population_2, type = "string", call = call)

  assert_vector(masked_markers, type = "character", n = 1, allow_null = TRUE, call = call)
  if (!is.null(masked_markers)) {
    assert_x_in_y(masked_markers, rownames(object), call = call)
  }
  assert_single_value(max_freq, type = "numeric", call = call)
  assert_within_limits(max_freq, limits = c(0, 1), call = call)
  assert_single_value(min_diff, type = "numeric", call = call)
  assert_within_limits(min_diff, limits = c(0, 1), call = call)
  assert_single_value(show_plot, type = "bool", call = call)

  assert_single_value(k, "integer", call = call)
  assert_within_limits(k, c(1, 6), call = call)
  assert_single_value(min_neighborhood_size, "integer", call = call)
  assert_within_limits(min_neighborhood_size, c(1, 50), call = call)
  assert_single_value(max_components_per_population, "integer", call = call)
  assert_within_limits(max_components_per_population, c(1, Inf), call = call)
  assert_single_value(neighborhoods_per_component, "integer", call = call)
  assert_within_limits(neighborhoods_per_component, c(10, 1e4), call = call)
  assert_single_value(verbose, type = "bool", call = call)

  return(invisible(NULL))
}

#' Utility function to score cell type identity of nodes using NNLS
#'
#' @param counts The count matrix for the nodes to be scored (rows are proteins,
#' columns are nodes/cells).
#' @param w The matrix of NMF weights for the two cell types (rows are proteins,
#' columns are cell types).
#'
#' @noRd
#'
.score_cell_types <- function(
  counts,
  w
) {
  shared_markers <- intersect(rownames(counts), rownames(w))
  x <- RcppML::project(
    A = counts[shared_markers, ] %>% as("dgCMatrix"),
    w = w[shared_markers, ],
    L1 = 0.2
  )

  x <- Matrix::t(x)

  x <- x + 1e-8
  x <- x %>%
    prop.table(margin = 1)
  rownames(x) <- colnames(counts)
  colnames(x) <- colnames(w)
  return(x)
}

#' Compute midpoint threshold from two k-means centers
#'
#' @param centers Numeric vector of length 2.
#'
#' @return Numeric midpoint between centers.
#'
#' @noRd
.kmeans_midpoint <- function(centers) {
  assert_vector(centers, type = "numeric", n = 2)
  mean(centers)
}

#' Perform spatial smoothing of a count matrix using the adjacency matrix of the graph
#'
#' @param A The adjacency matrix of the graph (sparse matrix of class `dgCMatrix`).
#' @param x The count matrix to be smoothed (sparse matrix of class `dgCMatrix`).
#' @param iter The number of iterations to perform the smoothing. Default is 5.
#'
#' @return A smoothed count matrix of the same dimensions as `x`.
#'
#' @export
#'
spatial_smoothing <- function(
  A,
  x,
  iter = 5L
) {
  assert_class(A, "dgCMatrix")
  assert_class(x, c("dgCMatrix", "matrix"))
  assert_singles_match(nrow(A), nrow(x))
  assert_single_value(iter, type = "integer")
  assert_within_limits(iter, limits = c(1, 50))

  # Assign equal weights to neighbors
  P <- A / Matrix::rowSums(A)

  # Assign center node the same weight as the sum of the neighbors
  Matrix::diag(P) <- 1

  # Renormalize to ensure rows sum to 1
  P <- P / Matrix::rowSums(P)

  # Apply smoothing
  for (i in seq_len(iter)) {
    x <- P %*% x
  }

  return(x)
}

#' Utility function to collapse a count matrix by a protein marker set
#'
#' In a PNA cell graph, each node is assigned a single protein label, meaning that the
#' node count matrix is binary. Either a node is positive for a protein (count = 1)
#' or it is not (count = 0). When collapsing the count matrix for a set of protein markers,
#' the count will always be either 1 or 0, so the returned count matrix is also binary.
#' This is essentially the same thing as treating a set of protein markers as a single marker.
#'
#' @param cg A `CellGraph` object containing the count matrix.
#' @param markers A character vector of protein names.
#'
#' @return A collapsed count matrix where counts for the specified markers are summed
#'
#' @noRd
#'
.create_collapsed_count_matrix <- function(
  cg,
  markers
) {
  x <- cg@counts[, markers, drop = FALSE] %>% Matrix::rowSums()
  m <- matrix(
    x,
    ncol = 1
  ) %>%
    as("dgCMatrix")
  rownames(m) <- rownames(cg@counts)
  return(m)
}

#' Utility function to subset a cell graph based on a boolean filter
#'
#' @param g An `tbl_graph` object representing a PNA cell graph.
#' @param bool_filter A boolean vector indicating which nodes to keep in the graph.
#' @param keep_largest_comp Whether to keep only the largest connected component after
#' subsetting. Default is TRUE.
#' @param min_comp_size Minimum size of connected components to keep if
#' `keep_largest_comp = FALSE`.
#' Note that the resulting graph can be disconnected.
#'
#' @return A filtered `tbl_graph` object
#'
#' @noRd
#'
.subset_graph <- function(
  g,
  bool_filter,
  keep_largest_comp = TRUE,
  min_comp_size = 10
) {
  if (sum(bool_filter) == length(g)) {
    return(g)
  }

  # Create membrane graph by subsetting on membrane nodes
  g_filtered <- g %>%
    filter(bool_filter)
  comps <- igraph::components(g_filtered)$membership
  comps_sizes <- table(comps)

  if (keep_largest_comp) {
    # Filter graph to include largest connected component
    g_filtered <- g_filtered %>%
      filter(comps == (which.max(comps_sizes) %>% names() %>% as.integer()))
  } else {
    comps_keep <- names(comps_sizes[comps_sizes >= min_comp_size]) %>% as.integer()
    g_filtered <- g_filtered %>%
      filter(comps %in% comps_keep)
  }

  return(g_filtered)
}

#' Utility function to classify nodes into two cell types using NNLS
#'
#' This function projects weights onto node neighborhood protein counts vectors to
#' quantify cell type identity. The classification is then derived by thresholding
#' the first population score at the midpoint between two k-means centers estimated
#' on that score distribution.
#'
#' @param counts_membrane A neighborhood count matrix for the nodes, with
#' dimensions mxp, where m is the number of nodes and p is the number of proteins.
#' @param A_k A reachability matrix for the k-hop neighborhood of the membrane nodes, with
#' dimensions mxm.
#' @param w A matrix of NMF weights for the two cell populations, with dimensions px2, where
#' p is the number of proteins and the columns represent the cell populations.
#' @param cell_names A character vector of length 2 containing the names of the two cell populations.
#' @param verbose Whether to print messages.
#'
#' @return A named vector of length m containing the cell type classification for each
#' node, where the names correspond to the node IDs and the values are either 1
#' (cell 1) or 2 (cell 2).
#'
#' @noRd
#'
.classify_nodes_nnls <- function(
  counts_membrane,
  A_k,
  A,
  w,
  cell_names,
  spatial_smoothing_iter = 5L,
  verbose = TRUE
) {
  # Classify cell1 and cell2 nodes using NNLS
  if (verbose) {
    cli::cli_alert_info(
      "Classifying {.str {cell_names}} nodes using NNLS"
    )
  }

  shared_markers <- intersect(colnames(counts_membrane), rownames(w))
  exp_counts <- A_k %*% counts_membrane[, shared_markers]
  # L1 normalization
  row_sums <- Matrix::rowSums(exp_counts)
  row_sums[row_sums == 0] <- 1
  exp_counts <- exp_counts / row_sums
  exp_counts <- as(exp_counts, "dgCMatrix") %>% Matrix::t()
  w <- w[shared_markers, , drop = FALSE]
  x <- RcppML::project(
    A = exp_counts,
    w = w,
    L1 = 0.2
  )
  x <- Matrix::t(x)

  if (spatial_smoothing_iter > 0) {
    if (verbose) {
      cli::cli_alert_info(
        "Performing spatial smoothing of projection scores with {spatial_smoothing_iter} iterations"
      )
    }
    x <- spatial_smoothing(A, x, iter = spatial_smoothing_iter)
    x <- as.matrix(x)
  }

  # Add a pseudocount to avoid getting non-finite values
  x <- x + 1e-8
  x <- x %>%
    prop.table(margin = 1)
  rownames(x) <- rownames(counts_membrane)
  colnames(x) <- colnames(w)

  # Slice scores into classes
  km <- kmeans(x[, cell_names[1]], centers = 2)
  th <- .kmeans_midpoint(km$centers)
  node_classification <- if_else(x[, cell_names[1]] >= th, 1, 2) %>%
    set_names(rownames(x))

  return(node_classification)
}

#' Utility function to fetch interface nodes between two cells in a graph
#'
#' An interface node is defined as a node that is connected to both cell 1 and cell 2
#' by a crossing edge. Once these interface nodes are identified, we can expand the
#' interface to include neighboring nodes.
#'
#' @param A The adjacency matrix of the graph, with dimensions xn, where n is the number of nodes.
#' @param c1_nodes A character vector of node IDs corresponding to cell 1.
#' @param c2_nodes A character vector of node IDs corresponding to cell 2.
#' @param k_interface_expansion An integer specifying how many hops to expand the interface nodes
#' to include neighboring nodes. Default is 1, meaning that we will include direct neighbors of
#' the interface nodes.
#'
#' @return A character vector of node IDs corresponding to the interface nodes.
#'
#' @noRd
#'
.fetch_interface_nodes <- function(
  A,
  c1_nodes,
  c2_nodes,
  k_interface_expansion = 1L
) {
  if (length(c1_nodes) == 0 || length(c2_nodes) == 0) {
    return(character(0))
  }

  # Slice out the adjacency matrix for crossing edges between cell1 and cell2
  A_crossing <- A[c1_nodes, c2_nodes, drop = FALSE]

  # Select the node IDs for nodes attached to a crossing edge
  interface_nodes <- c(
    rownames(A_crossing)[Matrix::rowSums(A_crossing) > 0],
    colnames(A_crossing)[Matrix::colSums(A_crossing) > 0]
  ) %>% unique()

  # Expand interface to include neighboring nodes.
  if (k_interface_expansion > 0) {
    A_k_interface <- expand_adjacency_matrix(A, k = k_interface_expansion)
    expanded_nodes <- colnames(A_k_interface)[Matrix::colSums(A_k_interface[interface_nodes, ]) > 0]
    interface_nodes <- unique(c(interface_nodes, expanded_nodes))
  }

  return(interface_nodes)
}

#' Compute distance from a set of seed nodes using fast breadth-first search.
#'
#' This function computes the shortest path distance from a set of seed nodes to all
#' other nodes in the graph using a breadth-first search (BFS) approach. It iteratively
#' expands the frontier of reached nodes until all reachable nodes have been assigned a
#' distance or the maximum number of iterations is reached.
#'
#' @param cg A `CellGraph` object.
#' @param seed_nodes A character vector of node names to use as seeds for distance calculation.
#' @param max_iter An integer specifying the maximum number of iterations (distance levels)
#' to compute. Default is 40.
#' @param verbose A logical value indicating whether to print progress messages during the
#' computation. Default is FALSE.
#'
#' @return A `CellGraph` object with an added `distance_from_seed` column in the node data,
#' indicating the shortest path distance from the nearest seed node.
#'
#' @examples
#' library(dplyr)
#' library(tidygraph)
#' cg <- ReadPNA_Seurat(minimal_pna_pxl_file(), verbose = FALSE) %>%
#'   LoadCellGraphs(cells = colnames(.)[4], add_layout = TRUE, verbose = FALSE) %>%
#'   CellGraphs() %>%
#'   .[[4]]
#'
#' # Compute distances from the seed set, here we just pick a random point
#' start_set <- cg@cellgraph %N>%
#'   pull(name) %>%
#'   head(1)
#' cg <- distance_from_node_set(cg, start_set)
#'
#' # Visualize the distance on the 3D layout
#' xyz <- cg@layout$wpmds_3d %>%
#'   mutate(d = cg@cellgraph %N>% pull(distance_from_seed))
#'
#' plotly::plot_ly(
#'   data = xyz,
#'   x = ~x, y = ~y, z = ~z,
#'   color = ~d,
#'   colors = c("lightgrey", "mistyrose", "red", "darkred", "black"),
#'   type = "scatter3d",
#'   mode = "markers",
#'   marker = list(size = 2)
#' )
#'
#' @export
#'
distance_from_node_set <- function(
  cg,
  seed_nodes,
  max_iter = 40L,
  verbose = FALSE
) {
  assert_class(cg, "CellGraph")
  assert_vector(seed_nodes, "character", n = 1)
  assert_single_value(max_iter, type = "integer")
  assert_within_limits(max_iter, limits = c(0, Inf))
  assert_single_value(verbose, type = "bool")

  nodes <- cg@cellgraph %N>% pull(name)
  n <- length(nodes)

  missing_seeds <- setdiff(seed_nodes, nodes)
  if (length(missing_seeds) > 0L) {
    cli::cli_abort(c(
      "x" = "All seed nodes must be present in the graph.",
      "i" = "The following seed nodes are not present in the graph: {missing_seeds}"
    ))
  }

  # Adjacency matrix — use column-major (dgCMatrix) and operate column-wise
  A <- igraph::as_adjacency_matrix(cg@cellgraph, sparse = TRUE)

  # Pre-allocate the distance vector (NA = unreached)
  d <- rep(NA_integer_, n)
  is_seed <- nodes %in% seed_nodes
  d[is_seed] <- 0L

  # Frontier = nodes added in the previous iteration (initially the seeds)
  frontier_idx <- which(is_seed)

  for (i in seq_len(max_iter)) {
    # Neighbors of the frontier: sum of frontier rows of A, restricted to unreached nodes
    f <- numeric(n)
    f[frontier_idx] <- 1
    neighbor_counts <- as.numeric(f %*% A)

    # New frontier = neighbors that have not been reached yet
    new_frontier_idx <- which(neighbor_counts > 0 & is.na(d))
    if (length(new_frontier_idx) == 0L) break

    d[new_frontier_idx] <- i
    frontier_idx <- new_frontier_idx

    if (isTRUE(verbose)) {
      cli::cli_inform("Iteration {i}: {length(new_frontier_idx)} new nodes reached.")
    }
  }

  # Add distance
  cg@cellgraph <- cg@cellgraph %N>%
    select(-any_of("distance_from_seed")) %>%
    left_join(tibble(name = nodes, distance_from_seed = d), by = "name")

  return(cg)
}

#' Compute partition counts
#'
#' @param cg A CellGraph object containing the cell graph and count data.
#' @param partition A character or factor vector indicating the partition of
#' nodes into groups. Either `partition` or `partition_column` must be provided.
#' @param partition_column A string indicating the name of the vertex attribute in
#' the cell graph that contains the partition information.
#'
#' @return A matrix of counts for each partition group, where rows correspond to
#' partition groups and columns correspond to features (e.g., proteins).
#'
#' @export
#'
partition_counts <- function(
  cg,
  partition = NULL,
  partition_column = NULL
) {
  assert_class(cg, "CellGraph")
  if (is.null(partition) && is.null(partition_column)) {
    cli::cli_abort("Either `partition` or `partition_column` must be provided.")
  }
  if (!is.null(partition) && !is.null(partition_column)) {
    cli::cli_abort("One of `partition` or `partition_column` must be provided, not both.")
  }
  assert_class(partition, c("character", "factor"), allow_null = TRUE)
  if (!is.null(partition)) {
    if (length(partition) != length(cg@cellgraph)) {
      cli::cli_abort(
        "Length of `partition` must match the number of nodes in the cell graph."
      )
    }
  }
  assert_single_value(partition_column, "string", allow_null = TRUE)

  if (!is.null(partition_column)) {
    partition <- try(
      cg@cellgraph %N>%
        pull(!!sym(partition_column)),
      silent = TRUE
    )
    if (inherits(partition, "try-error")) {
      cli::cli_abort(
        "Column '{partition_column}' not found in cell graph node attributes."
      )
    }
  }

  if (is.character(partition)) {
    lvls <- unique(partition)
    partition <- factor(partition, levels = lvls)
  }

  mm <- Matrix::sparse.model.matrix(~ 0 + partition)
  part_counts <- Matrix::t(mm) %*% cg@counts
  rownames(part_counts) <- levels(partition)
  return(part_counts)
}
