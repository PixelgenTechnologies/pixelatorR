# Convert objects to a [`CellGraphAssay`](CellGraphAssay-class.md)

Convert objects to a [`CellGraphAssay`](CellGraphAssay-class.md)

## Usage

``` r
as.CellGraphAssay(x, ...)

# S3 method for class 'Assay'
as.CellGraphAssay(
  x,
  cellgraphs = NULL,
  polarization = NULL,
  colocalization = NULL,
  fs_map = NULL,
  ...
)
```

## Arguments

- x:

  An object to convert to class
  [`CellGraphAssay`](CellGraphAssay-class.md)

- ...:

  Arguments passed to other methods

- cellgraphs:

  A list of [`CellGraph`](CellGraph-class.md) objects

- polarization:

  A `tbl_df` with polarization scores

- colocalization:

  A `tbl_df` with colocalization scores

- fs_map:

  A `tbl_df` with information on source pxl file paths, sample IDs, and
  component IDs

## Value

A `CellGraphAssay` object

## Examples

``` r
library(pixelatorR)
library(SeuratObject)
#> Loading required package: sp
#> 
#> Attaching package: 'SeuratObject'
#> The following objects are masked from 'package:base':
#> 
#>     intersect, t
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(tidygraph)
#> 
#> Attaching package: 'tidygraph'
#> The following object is masked from 'package:stats':
#> 
#>     filter

pxl_file <- minimal_mpx_pxl_file()
counts <- ReadMPX_counts(pxl_file)
#> ℹ Loading count data from C:/Users/max/AppData/Local/R/win-library/4.5/pixelatorR/extdata/five_cells/five_cells.pxl
edgelist <- ReadMPX_item(pxl_file, items = "edgelist")
#> ℹ Loading item(s) from: C:/Users/max/AppData/Local/R/win-library/4.5/pixelatorR/extdata/five_cells/five_cells.pxl
#> →   Loading edgelist data
#> ✔ Returning a 'tbl_df' object
components <- colnames(counts)
edgelist_split <-
  edgelist %>%
  select(upia, upib, component) %>%
  distinct() %>%
  group_by(component) %>%
  group_split() %>%
  setNames(nm = components)

# Convert data into a list of CellGraph objects
bipartite_graphs <- lapply(edgelist_split, function(x) {
  x <- x %>% as_tbl_graph(directed = FALSE)
  x <- x %>% mutate(node_type = case_when(name %in% edgelist$upia ~ "A", TRUE ~ "B"))
  attr(x, "type") <- "bipartite"
  CreateCellGraphObject(cellgraph = x)
})

# Create Assay
assay <- CreateAssayObject(counts = counts)

# Convert Assay to CellGraphAssay
cg_assay <- as.CellGraphAssay(assay, cellgraphs = bipartite_graphs)
cg_assay
#> CellGraphAssay data with 80 features for 5 cells
#> First 10 features:
#>  CD274, CD44, CD25, CD279, CD41, HLA-ABC, CD54, CD26, CD27, CD38 
#> Loaded CellGraph objects:
#>  5
```
