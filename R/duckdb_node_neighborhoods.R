#' Fetch node neighborhoods from pxl file(s) attached to a Seurat object
#'
#' @param object A Seurat object with a pna assay and attached pxl file(s)
#' @param components Optional character vector of component names to fetch
#' neighborhoods for. Defaults to all components in the object.
#' @param nodes_per_component Integer specifying how many start nodes to sample
#' per component. Defaults to 1000. If NULL, select all nodes in each component.
#' @param k Integer specifying how many hops to traverse in the neighborhood.
#' Defaults to 3.
#'
#' @examples
#' library(dplyr)
#' se <- ReadPNA_Seurat(minimal_pna_pxl_file())
#' node_nbs_matrix <- fetch_node_neighborhoods(
#'   se,
#'   nodes_per_component = 500,
#'   k = 2
#' )
#'
#' # Compare output with node neighborhoods from CellGraph
#' cg <- se %>%
#'   LoadCellGraphs(cells = colnames(.)[1], verbose = FALSE) %>%
#'   CellGraphs() %>%
#'   .[[1]]
#'
#' A <- igraph::as_adjacency_matrix(cg@cellgraph)
#' A2 <- expand_adjacency_matrix(A, k = 2)
#' Matrix::diag(A2) <- 1
#' node_nbs_matrix_cg <- Matrix::t(A2 %*% cg@counts)
#'
#' shared_nbs <- intersect(colnames(node_nbs_matrix), colnames(node_nbs_matrix_cg))
#' shared_markers <- intersect(rownames(node_nbs_matrix), rownames(node_nbs_matrix_cg))
#' node_nbs_matrix_cg <- node_nbs_matrix_cg[shared_markers, shared_nbs]
#' node_nbs_matrix <- node_nbs_matrix[shared_markers, shared_nbs]
#'
#' identical(node_nbs_matrix, node_nbs_matrix_cg)
#'
#' @return A sparse matrix where rows are proteins and columns are start nodes
#'
#' @export
#'
fetch_node_neighborhoods <- function(
  object,
  components = NULL,
  nodes_per_component = 1000,
  k = 3
) {
  assert_class(object, "Seurat")
  assay <- DefaultAssay(object)
  assert_pna_assay(object[[assay]])
  assert_vector(components, "character", n = 1, allow_null = TRUE)
  assert_x_in_y(components, colnames(object), allow_null = TRUE)
  assert_single_value(nodes_per_component, "integer", allow_null = TRUE)
  if (!is.null(nodes_per_component)) {
    assert_within_limits(nodes_per_component, c(1, 10000))
  }
  assert_single_value(k, "integer")
  assert_within_limits(k, c(1, 6))

  components <- components %||% colnames(object)

  el <- Edgelists(object, lazy = TRUE)

  node_nbs_matrix_list <- pbapply::pblapply(components, function(nm) {
    .fetch_node_neighborhoods_duckdb(el$src$con, nm, nodes_per_component, k) %>%
      .convert_nbs_table_to_matrix() %>%
      Matrix::t()
  })

  if (length(node_nbs_matrix_list) == 1) {
    node_nbs_matrix <- node_nbs_matrix_list[[1]]
  } else {
    node_nbs_matrix <- SeuratObject::RowMergeSparseMatrices(
      node_nbs_matrix_list[[1]],
      node_nbs_matrix_list[-1]
    )
  }

  return(node_nbs_matrix)
}

#' Fetch node neighborhoods from a pxl file
#'
#' @param con A DuckDB connection
#' @param component The component to fetch neighborhoods for
#' @param nodes_per_component The number of start nodes to sample per component
#' @param k The number of hops to traverse
#'
#' @noRd
#'
.fetch_node_neighborhoods_duckdb <- function(
  con,
  component,
  nodes_per_component,
  k
) {
  # If nodes_per_component is NULL, use 'ALL' for the SQL LIMIT clause
  limit_val <- if (is.null(nodes_per_component)) "ALL" else nodes_per_component

  node_nbs <- DBI::dbGetQuery(
    con,
    glue::glue("
        WITH RECURSIVE
        -- 1. Isolate the component's edges once at the very beginning
        comp_edges AS (
            SELECT umi1, umi2, marker_1, marker_2
            FROM combined_edgelist
            WHERE component = '{component}'
        ),
        -- 2. Create the proteins view from the isolated edges
        proteins_cte AS (
            SELECT umi1 AS name, marker_1 AS protein, 'umi1' AS original_role FROM comp_edges
            UNION
            SELECT umi2 AS name, marker_2 AS protein, 'umi2' AS original_role FROM comp_edges
        ),
        -- 3. Select seeds
        seeds AS (
            SELECT name AS start_node
            FROM proteins_cte
            GROUP BY name
            ORDER BY RANDOM()
            LIMIT {limit_val}
        ),
        -- 4. Define local adjacency
        local_adj AS (
            SELECT umi1 AS src, umi2 AS dst FROM comp_edges
            UNION ALL
            SELECT umi2 AS src, umi1 AS dst FROM comp_edges
        ),
        -- 5. Perform the k-hop traversal
        traversal AS (
            SELECT start_node, start_node AS current_node, 0 AS depth
            FROM seeds
            -- CRITICAL FIX: Use UNION instead of UNION ALL to prevent path explosion
            UNION
            SELECT t.start_node, l.dst, t.depth + 1
            FROM traversal t
            JOIN local_adj l ON t.current_node = l.src
            WHERE t.depth < {k}
        ),
        -- 6. Get distinct nodes reached per start_node
        distinct_nodes AS (
            SELECT DISTINCT start_node, current_node AS node_id
            FROM traversal
        )
        -- 7. Directly map reached nodes to their markers. Skip edge reconstruction!
        SELECT DISTINCT
            dn.start_node || '-' || p_start.original_role AS start_node,
            dn.node_id AS umi,
            p_reached.protein AS marker
        FROM distinct_nodes dn
        JOIN proteins_cte p_start ON dn.start_node = p_start.name
        JOIN proteins_cte p_reached ON dn.node_id = p_reached.name;
        ")
  )

  return(node_nbs)
}

#' Convert the node neighborhoods table to a sparse matrix
#'
#' @param node_nbs A data frame with columns start_node, umi, and marker
#'
#' @noRd
#'
.convert_nbs_table_to_matrix <- function(
  node_nbs
) {
  # SQL already returns distinct (start_node, umi, marker)
  # We just ensure it's clean before factoring
  node_nbs_counts <- node_nbs %>%
    dplyr::select(start_node, umi, marker) %>%
    dplyr::distinct()

  # Construct neighborhood count matrix
  sn_fct <- factor(node_nbs_counts$start_node, levels = unique(node_nbs_counts$start_node))
  protein_fct <- factor(node_nbs_counts$marker, levels = unique(node_nbs_counts$marker))

  node_nbs_count_matrix <- Matrix::sparseMatrix(
    i = as.integer(sn_fct),
    j = as.integer(protein_fct),
    x = 1, # Multiple UMIs with the same marker will sum together here
    dims = c(length(levels(sn_fct)), length(levels(protein_fct))),
    dimnames = list(levels(sn_fct), levels(protein_fct))
  )

  return(node_nbs_count_matrix)
}
