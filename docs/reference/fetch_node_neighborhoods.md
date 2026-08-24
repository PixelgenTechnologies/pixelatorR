# Fetch node neighborhoods from pxl file(s) attached to a Seurat object

Fetch node neighborhoods from pxl file(s) attached to a Seurat object

## Usage

``` r
fetch_node_neighborhoods(
  object,
  components = NULL,
  nodes_per_component = 1000,
  k = 3
)
```

## Arguments

- object:

  A Seurat object with a pna assay and attached pxl file(s)

- components:

  Optional character vector of component names to fetch neighborhoods
  for. Defaults to all components in the object.

- nodes_per_component:

  Integer specifying how many start nodes to sample per component.
  Defaults to 1000. If NULL, select all nodes in each component.

- k:

  Integer specifying how many hops to traverse in the neighborhood.
  Defaults to 3.

## Value

A sparse matrix where rows are proteins and columns are start nodes

## Examples

``` r
library(dplyr)
se <- ReadPNA_Seurat(minimal_pna_pxl_file())
#> ✔ Created a <Seurat> object with 5 cells and 158 targeted surface proteins
node_nbs_matrix <- fetch_node_neighborhoods(
  se,
  nodes_per_component = 500,
  k = 2
)

# Compare output with node neighborhoods from CellGraph
cg <- se %>%
  LoadCellGraphs(cells = colnames(.)[1], verbose = FALSE) %>%
  CellGraphs() %>%
  .[[1]]

A <- igraph::as_adjacency_matrix(cg@cellgraph)
A2 <- expand_adjacency_matrix(A, k = 2)
#> Warning: 'as(<dgCMatrix>, "ngCMatrix")' is deprecated.
#> Use 'as(., "nMatrix")' instead.
#> See help("Deprecated") and help("Matrix-deprecated").
Matrix::diag(A2) <- 1
node_nbs_matrix_cg <- Matrix::t(A2 %*% cg@counts)

shared_nbs <- intersect(colnames(node_nbs_matrix), colnames(node_nbs_matrix_cg))
shared_markers <- intersect(rownames(node_nbs_matrix), rownames(node_nbs_matrix_cg))
node_nbs_matrix_cg <- node_nbs_matrix_cg[shared_markers, shared_nbs]
node_nbs_matrix <- node_nbs_matrix[shared_markers, shared_nbs]

identical(node_nbs_matrix, node_nbs_matrix_cg)
#> [1] TRUE
```
