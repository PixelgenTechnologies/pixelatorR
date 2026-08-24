# Create a CellGraph object

Create a CellGraph object

## Usage

``` r
CreateCellGraphObject(cellgraph, counts = NULL, layout = NULL, verbose = FALSE)
```

## Arguments

- cellgraph:

  A `tbl_graph` object representing a PNA single-cell graph

- counts:

  A `dgCMatrix` with marker counts

- layout:

  A named `list` of `tbl_df` objects with cell layouts

- verbose:

  Print messages

## Value

A `CellGraph` object

## Examples

``` r

library(pixelatorR)
library(dplyr)
library(tidygraph)

# Open a database connection (PXL file)
db <- PixelDB$new(minimal_pna_pxl_file())

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
node_names <- component_graph %>% pull(name)
# Ensure that the counts matrix rows match the graph node names
counts <- counts[node_names, ]

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
# Ensure that the layout table rows match the graph node names
layout <- layout[match(node_names, layout$name), ] %>%
  select(-name)

# Create a CellGraph object with graph, counts and layout
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
