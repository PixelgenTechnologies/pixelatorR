# Create a CellGraphAssay object

Create a [`CellGraphAssay`](CellGraphAssay-class.md) object from a count
matrix. The expected format of the input matrix is features x cells.

## Usage

``` r
CreateCellGraphAssay(
  counts,
  cellgraphs,
  polarization = NULL,
  colocalization = NULL,
  fs_map = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- counts:

  Unnormalized data (raw counts)

- cellgraphs:

  A named list of [`CellGraph`](CellGraph-class.md) objects

- polarization:

  A `tbl_df` with polarization scores

- colocalization:

  A `tbl_df` with colocalization scores

- fs_map:

  A `tbl_df` with information on source pxl file paths, sample IDs, and
  component IDs

- verbose:

  Print messages

- ...:

  Additional arguments passed to
  [`CreateAssayObject`](https://satijalab.github.io/seurat-object/reference/CreateAssayObject.html)

## Value

A `CellGraphAssay` object

## Examples

``` r
library(pixelatorR)
library(dplyr)
library(tidygraph)

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

# Create CellGraphAssay
cg_assay <- CreateCellGraphAssay(counts = counts, cellgraphs = bipartite_graphs)
cg_assay
#> CellGraphAssay data with 80 features for 5 cells
#> First 10 features:
#>  CD274, CD44, CD25, CD279, CD41, HLA-ABC, CD54, CD26, CD27, CD38 
#> Loaded CellGraph objects:
#>  5
```
