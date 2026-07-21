# PNAAssay5 Methods

Methods for [`PNAAssay5`](PNAAssay5-class.md) objects for generics
defined in other packages

Join Layers Together

## Usage

``` r
# S3 method for class 'PNAAssay5'
RenameCells(object, new.names = NULL, ...)

# S4 method for class 'PNAAssay5'
show(object)

# S3 method for class 'PNAAssay5'
subset(x, features = NULL, cells = NULL, ...)

# S3 method for class 'PNAAssay5'
merge(
  x = NULL,
  y = NULL,
  merge.data = TRUE,
  add.cell.ids = NULL,
  collapse = TRUE,
  ...
)

# S3 method for class 'PNAAssay5'
JoinLayers(object, layers = NULL, new = NULL, ...)
```

## Arguments

- object:

  A `PNAAssay5` object.

- new.names:

  A character vector with new cell IDs. The length of the vector must be
  equal to the number of cells in the object and the names must be
  unique.

- ...:

  Additional arguments passed to other methods

- x:

  A [`PNAAssay5`](PNAAssay5-class.md) object

- features:

  Feature names

- cells:

  Cell names

- y:

  A [`PNAAssay5`](PNAAssay5-class.md) object or a list of
  [`PNAAssay5`](PNAAssay5-class.md) objects

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

A `PNAAssay5` object with layers joined

## Functions

- `RenameCells(PNAAssay5)`: Rename cell IDs of a `PNAAssay5` object

- `show(PNAAssay5)`: Show method for `PNAAssay5` objects

- `subset(PNAAssay5)`: Subset a `PNAAssay5` object

- `merge(PNAAssay5)`: Merge two or more `PNAAssay5` objects together

- `JoinLayers(PNAAssay5)`: Join layers

## Examples

``` r
library(SeuratObject)
options(Seurat.object.assay.version = "v5")
# Load example data as a Seurat object
pxl_file <- minimal_pna_pxl_file()
seur_obj <- ReadPNA_Seurat(pxl_file)

# Merge Seurat objects
seur_obj_merged <- merge(seur_obj, seur_obj)
#> Warning: Some cell names are duplicated across objects provided. Renaming to enforce unique cell names.
pna_assay <- seur_obj_merged[["PNA"]]

# The PNAAssay5 now has two count matrices
pna_assay
#> PNAAssay (v5) data with 158 features for 10 cells
#> Top 10 variable features:
#>  HLA-ABC, B2M, CD11b, CD11c, CD18, CD82, CD8, TCRab, HLA-DR, CD45 
#> Layers:
#>  counts.1, counts.2 
#> Loaded CellGraph objects:
#>  0

# Join layers
pna_assay <- JoinLayers(pna_assay)

# Now the PNAAssay5 has a single, merged count matrix
pna_assay
#> PNAAssay (v5) data with 158 features for 10 cells
#> Top 10 variable features:
#>  HLA-ABC, B2M, CD11b, CD11c, CD18, CD82, CD8, TCRab, HLA-DR, CD45 
#> Layers:
#>  counts 
#> Loaded CellGraph objects:
#>  0

# JoinLayers now also works on the Seurat object directly
seur_obj_merged <- JoinLayers(seur_obj_merged)
seur_obj_merged[["PNA"]]
#> PNAAssay (v5) data with 158 features for 10 cells
#> Top 10 variable features:
#>  HLA-ABC, B2M, CD11b, CD11c, CD18, CD82, CD8, TCRab, HLA-DR, CD45 
#> Layers:
#>  counts 
#> Loaded CellGraph objects:
#>  0
```
