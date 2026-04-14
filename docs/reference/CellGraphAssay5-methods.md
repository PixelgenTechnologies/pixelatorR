# CellGraphAssay5 Methods

Methods for [`CellGraphAssay5`](CellGraphAssay5-class.md) objects for
generics defined in other packages

Join Layers Together

## Usage

``` r
# S3 method for class 'CellGraphAssay5'
RenameCells(object, new.names = NULL, ...)

# S4 method for class 'CellGraphAssay5'
show(object)

# S3 method for class 'CellGraphAssay5'
subset(x, features = NULL, cells = NULL, ...)

# S3 method for class 'CellGraphAssay5'
merge(
  x = NULL,
  y = NULL,
  merge.data = TRUE,
  add.cell.ids = NULL,
  collapse = TRUE,
  ...
)

# S3 method for class 'CellGraphAssay5'
JoinLayers(object, layers = NULL, new = NULL, ...)
```

## Arguments

- object:

  A `CellGraphAssay5` object.

- new.names:

  A character vector with new cell IDs. The length of the vector must be
  equal to the number of cells in the object and the names must be
  unique.

- ...:

  Additional arguments passed to other methods

- x:

  A [`CellGraphAssay5`](CellGraphAssay5-class.md) object

- features:

  Feature names

- cells:

  Cell names

- y:

  A [`CellGraphAssay5`](CellGraphAssay5-class.md) object or a list of
  [`CellGraphAssay5`](CellGraphAssay5-class.md) objects

- merge.data:

  Merge the data slots instead of just merging the counts (which
  requires renormalization); this is recommended if the same
  normalization approach was applied to all objects

- add.cell.ids:

  A character vector with sample names

- collapse:

  If TRUE, merge layers of the same name together

- layers:

  A character vector of layer names to join.

- new:

  Name of new layers

## Value

A `CellGraphAssay5` object with layers joined

## Functions

- `RenameCells(CellGraphAssay5)`: Rename cell IDs of a `CellGraphAssay5`
  object

- `show(CellGraphAssay5)`: Show method for `CellGraphAssay5` objects

- `subset(CellGraphAssay5)`: Subset a `CellGraphAssay5` object

- `merge(CellGraphAssay5)`: Merge two or more `CellGraphAssay5` objects
  together

- `JoinLayers(CellGraphAssay5)`: Join layers

## Examples

``` r
library(SeuratObject)
options(Seurat.object.assay.version = "v5")

# Load example data as a Seurat object
pxl_file <- minimal_mpx_pxl_file()
seur_obj <- ReadMPX_Seurat(pxl_file)
#> ✔ Created a 'Seurat' object with 5 cells and 80 targeted surface proteins
#> ! Failed to remove temporary file C:/Users/max/AppData/Local/Temp/RtmpmOhqBt/file62e864d31e21.h5ad

# Merge Seurat objects
seur_obj_merged <- merge(seur_obj, seur_obj)
#> Warning: Some cell names are duplicated across objects provided. Renaming to enforce unique cell names.
cg_assay <- seur_obj_merged[["mpxCells"]]

# The CellGraphAssay5 now has two count matrices
cg_assay
#> CellGraphAssay (v5) data with 80 features for 10 cells
#> Top 10 variable features:
#>  CD274, CD44, CD25, CD279, CD41, HLA-ABC, CD54, CD26, CD27, CD38 
#> Layers:
#>  counts.1, counts.2 
#> Loaded CellGraph objects:
#>  0

# Join layers
cg_assay <- JoinLayers(cg_assay)

# Now the CellGraphAssay5 has a single, merged count matrix
cg_assay
#> CellGraphAssay (v5) data with 80 features for 10 cells
#> Top 10 variable features:
#>  CD274, CD44, CD25, CD279, CD41, HLA-ABC, CD54, CD26, CD27, CD38 
#> Layers:
#>  counts 
#> Loaded CellGraph objects:
#>  0

# JoinLayers now also works on the Seurat object directly
seur_obj_merged <- JoinLayers(seur_obj_merged)
seur_obj_merged[["mpxCells"]]
#> CellGraphAssay (v5) data with 80 features for 10 cells
#> Top 10 variable features:
#>  CD274, CD44, CD25, CD279, CD41, HLA-ABC, CD54, CD26, CD27, CD38 
#> Layers:
#>  counts 
#> Loaded CellGraph objects:
#>  0
```
