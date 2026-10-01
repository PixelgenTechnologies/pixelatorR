# Show method for `CellGraph` object

Show method for `CellGraph` object

## Usage

``` r
# S4 method for class 'CellGraph'
show(object)
```

## Arguments

- object:

  A `CellGraph` object

## Examples

``` r
library(pixelatorR)
library(dplyr)
library(tidygraph)

edge_list <-
  ReadMPX_item(
    minimal_mpx_pxl_file(),
    items = "edgelist"
  )
bipart_graph <-
  edge_list %>%
  select(upia, upib, marker) %>%
  distinct() %>%
  as_tbl_graph(directed = FALSE) %>%
  mutate(node_type = case_when(name %in% edge_list$upia ~ "A", TRUE ~ "B"))
attr(bipart_graph, "type") <- "bipartite"

cg <- CreateCellGraphObject(cellgraph = bipart_graph)

# Show method
cg
#> A CellGraph object containing a bipartite graph with 16800 nodes and 68255 edges
```
