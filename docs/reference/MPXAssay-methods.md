# MPXAssay Methods

Methods for [`MPXAssay`](MPXAssay-class.md) objects for generics defined
in other packages

## Usage

``` r
# S3 method for class 'MPXAssay'
RenameCells(object, new.names = NULL, ...)

# S4 method for class 'MPXAssay'
show(object)

# S3 method for class 'MPXAssay'
subset(x, features = NULL, cells = NULL, ...)

# S3 method for class 'MPXAssay'
merge(
  x = NULL,
  y = NULL,
  merge.data = TRUE,
  add.cell.ids = NULL,
  collapse = TRUE,
  ...
)
```

## Arguments

- object:

  A `CellGraphAssay` or a `CellGraphAssay5` object

- new.names:

  A character vector with new cell IDs. The length of the vector must be
  equal to the number of cells in the object and the names must be
  unique.

- ...:

  Arguments passed to other methods

- x:

  A [`MPXAssay`](MPXAssay-class.md) object

- features:

  Feature names

- cells:

  Cell names

- y:

  A [`MPXAssay`](MPXAssay-class.md) object or a list of
  [`MPXAssay`](MPXAssay-class.md) objects

- merge.data:

  Merge the data slots instead of just merging the counts (which
  requires renormalization); this is recommended if the same
  normalization approach was applied to all objects

- add.cell.ids:

  A character vector with sample names

- collapse:

  If TRUE, merge layers of the same name together

## Value

A `MPXAssay` object

## Functions

- `RenameCells(MPXAssay)`: Rename cell IDs of a `CellGraphAssay` or
  `CellGraphAssay5` object

- `show(MPXAssay)`: Show method for a `CellGraphAssay` or a
  `CellGraphAssay5` object

- `subset(MPXAssay)`: Subset a `CellGraphAssay` or a `CellGraphAssay5`
  object

- `merge(MPXAssay)`: Merge two or more `CellGraphAssay` or
  `CellGraphAssay5` objects together

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

# Show method
cg_assay
#> CellGraphAssay data with 80 features for 5 cells
#> First 10 features:
#>  CD274, CD44, CD25, CD279, CD41, HLA-ABC, CD54, CD26, CD27, CD38 
#> Loaded CellGraph objects:
#>  5

library(pixelatorR)
library(dplyr)
options(Seurat.object.assay.version = "v3")

pxl_file <- minimal_mpx_pxl_file()
seur <- ReadMPX_Seurat(pxl_file)
#> ✔ Created a 'Seurat' object with 5 cells and 80 targeted surface proteins
#> ! Failed to remove temporary file C:/Users/max/AppData/Local/Temp/RtmpmOhqBt/file62e851be3376.h5ad
seur <- LoadCellGraphs(seur)
#> →    Loading CellGraphs for 5 cells from sample 1
#> ✔ Successfully loaded 5 CellGraph object(s).
cg_assay <- seur[["mpxCells"]]

# Subset CellGraphAssay(5)
# ---------------------------------
cg_assay_subset <- subset(cg_assay, cells = colnames(cg_assay)[1:3])

# Subset Seurat object containing a CellGraphAssay(5)
# --------------------------------
seur_subset <- subset(seur, cells = colnames(seur)[1:3])

# Merge multiple CellGraphAssay(5) objects
# ---------------------------------

# Merge 3 CellGraphAssay(5) objects
cg_assay_merged <- merge(cg_assay,
  y = list(cg_assay, cg_assay),
  add.cell.ids = c("A", "B", "C")
)
```
