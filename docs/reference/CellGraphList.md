# The CellGraphList class

A named list of [`CellGraph`](CellGraph-class.md) objects. Unloaded
graphs may be stored as `NULL`. See
[`CellGraphList-methods`](CellGraphList-methods.md) for subsetting,
replacement, and concatenation.

## Usage

``` r
CreateCellGraphList(cellgraphs = list())
```

## Arguments

- cellgraphs:

  A named list of [`CellGraph`](CellGraph-class.md) objects. Unloaded
  graphs may be represented as `NULL`.

## Value

A `CellGraphList` object

## See also

[`CellGraphList-methods`](CellGraphList-methods.md)

## Examples

``` r
library(pixelatorR)
library(tidygraph)
library(dplyr)

# Build a small dummy cell graph
edges <- tibble(from = c("a", "b"), to = c("b", "c"))
g <- as_tbl_graph(edges, directed = FALSE) %N>%
  mutate(node_type = c("umi1", "umi2", "umi1"))
attr(g, "type") <- "bipartite"
cg <- CreateCellGraphObject(cellgraph = g)

# Repeat the CellGraph in a named list and convert to a CellGraphList
cgl <- CreateCellGraphList(list(cell_1 = cg, cell_2 = cg))
cgl
#> A CellGraphList with 2 loaded CellGraph object(s) out of 2 
#> Names: cell_1, cell_2  
```
