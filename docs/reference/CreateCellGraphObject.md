# Create a CellGraph object

Create a CellGraph object

## Usage

``` r
CreateCellGraphObject(
  cellgraph,
  counts = NULL,
  layout = NULL,
  layers = NULL,
  meta.data = NULL,
  reductions = NULL,
  verbose = FALSE
)
```

## Arguments

- cellgraph:

  A `tbl_graph` object representing a PNA single-cell graph

- counts:

  A `dgCMatrix` with marker counts. Rows are matched to graph node names
  (order does not need to match).

- layout:

  A named `list` of `data.frame` objects with cell layouts. Nodes are
  identified by row names or by a `name` column; otherwise the row order
  is assumed to follow the graph. MPX bipartite layouts may use
  unsuffixed names while graph nodes keep `-A`/`-B`; those names are
  matched after stripping the suffix, as in
  [`LoadCellGraphs`](LoadCellGraphs.md). Stored layouts keep graph node
  order and do not copy node IDs as row names.

- layers:

  A named `list` of additional numeric node matrices (nodes x features).
  `"counts"` is reserved.

- meta.data:

  A node-level `data.frame` or `tbl_df`. Either row names or a `name`
  column must identify nodes.

- reductions:

  A named `list` of [`NodeDimReduc`](NodeDimReduc-class.md) objects

- verbose:

  Print messages

## Value

A `CellGraph` object

## Details

Node-level variable names must not clash between the graph node table,
`meta.data`, reduction embeddings, and matrix features. Count and layer
matrices may share feature names because methods such as
[`FetchData`](https://satijalab.github.io/seurat-object/reference/FetchData.html)
select a specific layer.

## Examples

``` r
library(pixelatorR)
library(dplyr)
library(tidygraph)

# Open a database connection (PXL file)
db <- PixelDB$new(minimal_pna_pxl_file())
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpmC3mql/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.

# Select a component ID and load the edgelist
sel_comp <- db$cell_meta() %>%
  rownames() %>%
  head(1)
component_edgelist <- db$components_edgelist(
  components = sel_comp,
  umi_data_type = "suffixed_string"
) %>%
  select(umi1, umi2)

# Define node types for the bipartite graph
umi_node_type <- bind_rows(
  component_edgelist %>% select(name = umi1) %>% mutate(node_type = "umi1"),
  component_edgelist %>% select(name = umi2) %>% mutate(node_type = "umi2")
) %>%
  distinct()

# Create a bipartite graph from the edgelist and add node types
component_graph <- as_tbl_graph(component_edgelist, directed = FALSE) %N>%
  left_join(umi_node_type, by = "name")

# Set the graph type attribute to "bipartite"
attr(component_graph, "type") <- "bipartite"

# Create a CellGraph object with just the graph
cg <- CreateCellGraphObject(cellgraph = component_graph)
cg
#> A CellGraph object containing a bipartite graph with 43543 nodes and 97014 edges

# Load cell count matrix
counts <- db$components_marker_counts(
  components = sel_comp, as_sparse = TRUE
)[[1]]

# Create a CellGraph object with graph and counts
cg <- CreateCellGraphObject(cellgraph = component_graph, counts = counts)
cg
#> A CellGraph object containing a bipartite graph with 43543 nodes and 97014 edges
#> Number of markers:  149 

# Create a CellGraph object with counts and layout
layout <- db$components_layout(
  components = sel_comp
)[[1]]
#> ℹ Fetching 1 component layouts...

# Layouts with a name column or node row names are matched automatically
cg <- CreateCellGraphObject(
  cellgraph = component_graph,
  counts = counts,
  layout = list(wpmds_3d = layout)
)
cg
#> A CellGraph object containing a bipartite graph with 43543 nodes and 97014 edges
#> Number of markers:  149 
#> Layouts: wpmds_3d 
```
