# Create a CellGraph object

Create a CellGraph object

## Usage

``` r
CreateCellGraphObject(cellgraph, counts = NULL, layout = NULL, verbose = FALSE)
```

## Arguments

- cellgraph:

  A `tbl_graph` object representing an mpx single-cell graph

- counts:

  A `dgCMatrix` with marker counts

- layout:

  A `tbl_df` object with cell layout(s)

- verbose:

  Print messages

## Value

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
#> ℹ Loading item(s) from: C:/Users/max/AppData/Local/Temp/Rtmp2f7Zq5/temp_libpath7aa0ad0713/pixelatorR/extdata/five_cells/five_cells.pxl
#> →   Loading edgelist data
#> ✔ Returning a 'tbl_df' object
bipart_graph <-
  edge_list %>%
  select(upia, upib, marker) %>%
  distinct() %>%
  as_tbl_graph(directed = FALSE) %>%
  mutate(node_type = case_when(name %in% edge_list$upia ~ "A", TRUE ~ "B"))
attr(bipart_graph, "type") <- "bipartite"

cg <- CreateCellGraphObject(cellgraph = bipart_graph)
cg
#> A CellGraph object containing a bipartite graph with 16800 nodes and 68255 edges
```
