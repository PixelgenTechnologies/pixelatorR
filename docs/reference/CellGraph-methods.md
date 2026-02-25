# CellGraph Methods

Methods for [`CellGraph`](CellGraph-class.md) objects for generics
defined in other packages

## Usage

``` r
# S4 method for class 'CellGraph'
show(object)

# S3 method for class 'CellGraph'
subset(x, nodes, ...)
```

## Arguments

- object:

  A [`CellGraph`](CellGraph-class.md) object

- x:

  A [`CellGraph`](CellGraph-class.md) object

- nodes:

  A character vector of node names

- ...:

  Currently not used

## Value

A `CellGraph` object containing only the specified nodes.

## Functions

- `show(CellGraph)`: Show a `CellGraph` object

- `subset(CellGraph)`: Subset a `CellGraph` object

## Examples

``` r
library(pixelatorR)
se <- ReadPNA_Seurat(minimal_pna_pxl_file(), verbose = FALSE)
se <- LoadCellGraphs(se, cells = colnames(se)[1], verbose = FALSE)
cg <- CellGraphs(se)[[1]]

# Show method
cg
#> A CellGraph object containing a bipartite graph with 43543 nodes and 97014 edges
#> Number of markers:  149 

# Subset
cg_small <- subset(cg, nodes = rownames(cg@counts)[1:100])
cg_small
#> A CellGraph object containing a bipartite graph with 100 nodes and 53 edges
#> Number of markers:  149 
```
