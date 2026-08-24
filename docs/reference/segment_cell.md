# cell:cell conjugate segmentation

This function performs segmentation of a cell:cell conjugate graph into
its respective cell types and the interface between them (optional).

## Usage

``` r
segment_cell(
  cg,
  w,
  k = 2L,
  detect_interface = TRUE,
  k_interface_expansion = 0L,
  keep_largest_comp = TRUE,
  min_comp_size = 10,
  spatial_smoothing_iter = 5L,
  verbose = TRUE
)
```

## Arguments

- cg:

  A `CellGraph` object

- w:

  A matrix of NMF weights with rows corresponding to proteins and
  columns corresponding to the two interacting cell types.

- k:

  The number of hops to use when computing local node abundance profiles
  used to compute the NNLS projection scores.

- detect_interface:

  Whether to identify the interface between the two cell types.

- k_interface_expansion:

  The number of hops to expand the interface between the two cell types.
  The default (0) classifies nodes that are directly connected to both
  cell types as interface. Values greater than 0 can be set to expand
  the interface further by including nodes that are within k hops of the
  initial interface nodes.

- keep_largest_comp:

  Whether to keep only the largest connected component in each segmented
  cell-type graph. If `FALSE`, all connected components with at least
  `min_comp_size` nodes are kept.

- min_comp_size:

  The minimum size of connected components to keep when
  `keep_largest_comp` is `FALSE`.

- spatial_smoothing_iter:

  The number of iterations to perform for spatial smoothing of node
  neighborhoods scores. This helps smoothen out local noise in the graph
  and improve classification.

- verbose:

  Print updates during the segmentation process.

## Value

A `CellGraph` object containing the compartment labels from the
segmentation as a node attribute.

## Details

Requirements:

- A `CellGraph` object containing the graph representation of the
  cell:cell conjugate. Only two cell types should be present in the
  graph.

- A matrix of NMF weights (`w`) for the two cell types, which can be
  derived using the `cc_protein_weights` function.

The segmentation is performed in several steps:

1.  Classify nodes into "cell1" and "cell2" using Non-Negative Least
    Squares (NNLS). This method is described in detail below. After
    classification, the Largest Connected Component (LCC) is kept for
    each cell type to remove small disconnected subgraphs.

2.  (Optional) Identify the interface nodes between the two cell types.
    These are defined as nodes that are directly connected to both cell
    types. Optionally, the interface can be expanded by including
    neighboring nodes within `k_interface_expansion` hops. By setting
    `detect_interface` to `FALSE`, the interface detection step is
    skipped.

Nodes that cannot be classified with confidence are labeled as "other".

## NNLS method for classifying membrane nodes

To segment cell:cell conjugates, we first need to fit an NMF model to
the abundance matrix (or local node neighborhood abundance profiles) of
the two cell populations using the `cc_protein_weights` function. This
gives us the protein weights for each cell type (`w`) that we can use
for classification. As a rule of thumb, a good NMF model should have
high weights for proteins that are mutually exclusive between the two
cell types. For instance, if we are segmenting a B cell:T cell
conjugate, we would expect the NMF model to give high weights to CD10,
CD20 and CD21 for the B cell identity and CD3e, CD2 and CD5 for the T
cell identity. For the membrane segmentation step, we create a
neighborhood profile for each node by summing the abundance of its
neighbors in the graph. We then project the NMF weights (`w`) onto these
neighborhood profiles using NNLS to get a score for each node and cell
type.

## Examples

``` r
library(dplyr)
se <- ReadPNA_Seurat(minimal_pna_pxl_file()) %>%
  LoadCellGraphs(cells = colnames(.)[2], verbose = FALSE)
cg <- CellGraphs(se)[[2]]

# Get cell type weights
se$cell_type <- c("Mono", "pDC", "CD4T", "CD4T", "CD4T")
w <- cc_protein_weights(
  se,
  group_by = "cell_type",
  population_1 = "Mono", population_2 = "CD4T",
  show_plot = FALSE
)
#> ! Population "Mono" has less than 20 cells, which may lead to unstable results.
#> ! Population "CD4T" has less than 20 cells, which may lead to unstable results.

cg <- segment_cell(
  cg,
  w = w
)
#> ℹ Classifying "Mono" and "CD4T" nodes using NNLS
#> ℹ Performing spatial smoothing of projection scores with 5 iterations
#> ℹ Defining interface nodes between "Mono" and "CD4T"
#> ℹ Mapping nodes to compartments
```
