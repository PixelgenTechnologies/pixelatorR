# Convert objects to a [`CellGraphAssay5`](CellGraphAssay5-class.md)

Convert objects to a [`CellGraphAssay5`](CellGraphAssay5-class.md)

## Usage

``` r
as.CellGraphAssay5(x, ...)

# S3 method for class 'Assay5'
as.CellGraphAssay5(
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
  [`CellGraphAssay5`](CellGraphAssay5-class.md)

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

A `CellGraphAssay5` object

## Examples

``` r
library(pixelatorR)
library(SeuratObject)
library(dplyr)
library(tidygraph)

pxl_file <- minimal_mpx_pxl_file()
counts <- ReadMPX_counts(pxl_file)
#> ℹ Loading count data from C:/Users/max/AppData/Local/Temp/Rtmp2f7Zq5/temp_libpath7aa0ad0713/pixelatorR/extdata/five_cells/five_cells.pxl
edgelist <- ReadMPX_item(pxl_file, items = "edgelist")
#> ℹ Loading item(s) from: C:/Users/max/AppData/Local/Temp/Rtmp2f7Zq5/temp_libpath7aa0ad0713/pixelatorR/extdata/five_cells/five_cells.pxl
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

# Create Assay5
assay5 <- CreateAssay5Object(counts = counts)
#> Warning: Data is of class matrix. Coercing to dgCMatrix.

# Convert Assay5 to CellGraphAssay5
cg_assay <- as.CellGraphAssay5(assay5, cellgraphs = bipartite_graphs)
cg_assay
#> CellGraphAssay (v5) data with 80 features for 5 cells
#> First 10 features:
#>  CD274, CD44, CD25, CD279, CD41, HLA-ABC, CD54, CD26, CD27, CD38 
#> Layers:
#>  counts 
#> Loaded CellGraph objects:
#>  5
```
