# PNAAssay Methods

Methods for [`PNAAssay`](PNAAssay-class.md) objects for generics defined
in other packages

## Usage

``` r
# S3 method for class 'PNAAssay'
RenameCells(object, new.names = NULL, ...)

# S4 method for class 'PNAAssay'
show(object)

# S3 method for class 'PNAAssay'
subset(x, features = NULL, cells = NULL, ...)

# S3 method for class 'PNAAssay'
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

  A `PNAAssay` or a `PNAAssay5` object

- new.names:

  A character vector with new cell IDs. The length of the vector must be
  equal to the number of cells in the object and the names must be
  unique.

- ...:

  Arguments passed to other methods

- x:

  A [`PNAAssay`](PNAAssay-class.md) object

- features:

  Feature names

- cells:

  Cell names

- y:

  A [`PNAAssay`](PNAAssay-class.md) object or a list of
  [`PNAAssay`](PNAAssay-class.md) objects

- merge.data:

  Merge the data slots instead of just merging the counts (which
  requires renormalization); this is recommended if the same
  normalization approach was applied to all objects

- add.cell.ids:

  A character vector with sample names

- collapse:

  If TRUE, merge layers of the same name together

## Value

A `PNAAssay` object

## Functions

- `RenameCells(PNAAssay)`: Rename cell IDs of a `PNAAssay` or
  `PNAAssay5` object

- `show(PNAAssay)`: Show method for `PNAAssay` objects

- `subset(PNAAssay)`: Subset a `PNAAssay` or a `PNAAssay5` object

- `merge(PNAAssay)`: Merge two or more `PNAAssay` or `PNAAssay5` objects
  together

## Examples

``` r
library(pixelatorR)
library(SeuratObject)

pxl_file <- minimal_pna_pxl_file()
seur_obj <- ReadPNA_Seurat(pxl_file)
pna_assay <- seur_obj[["PNA"]]
pna_assay <- RenameCells(pna_assay, new.names = paste0(colnames(pna_assay), "-new"))
colnames(pna_assay)
#> [1] "0a45497c6bfbfb22-new" "2708240b908e2eba-new" "c3c393e9a17c1981-new"
#> [4] "d4074c845bb62800-new" "efe0ed189cb499fc-new"

library(pixelatorR)

pxl_file <- minimal_pna_pxl_file()
seur_obj <- ReadPNA_Seurat(pxl_file)
seur_obj[["PNA"]]
#> PNAAssay data with 158 features for 5 cells
#> Top 10 variable features:
#>  HLA-ABC, B2M, CD11b, CD11c, CD18, CD82, CD8, TCRab, HLA-DR, CD45 
#> Loaded CellGraph objects:
#>  0

library(pixelatorR)

pxl_file <- minimal_pna_pxl_file()
seur_obj <- ReadPNA_Seurat(pxl_file)
pna_assay <- seur_obj[["PNA"]]

pna_assay <- subset(pna_assay, cells = colnames(pna_assay)[1:2])
pna_assay
#> PNAAssay data with 158 features for 2 cells
#> Top 10 variable features:
#>  HLA-ABC, B2M, CD11b, CD11c, CD18, CD82, CD8, TCRab, HLA-DR, CD45 
#> Loaded CellGraph objects:
#>  0

library(pixelatorR)
library(dplyr)

pxl_file <- minimal_pna_pxl_file()
seur_obj <- ReadPNA_Seurat(pxl_file)
pna_assay <- seur_obj[["PNA"]]

# Merge two data sets
pna_assay_merged <-
  merge(pna_assay, pna_assay, add.cell.ids = c("A", "B"))
pna_assay_merged
#> PNAAssay data with 158 features for 10 cells
#> First 10 features:
#>  HLA-ABC, B2M, CD11b, CD11c, CD18, CD82, CD8, TCRab, HLA-DR, CD45 
#> Loaded CellGraph objects:
#>  0

# Merge multiple data sets
pna_assay_list <- list(pna_assay, pna_assay, pna_assay)
pna_assay_merged <-
  merge(pna_assay_list[[1]], pna_assay_list[-1], add.cell.ids = c("A", "B", "C"))
pna_assay_merged
#> PNAAssay data with 158 features for 15 cells
#> First 10 features:
#>  HLA-ABC, B2M, CD11b, CD11c, CD18, CD82, CD8, TCRab, HLA-DR, CD45 
#> Loaded CellGraph objects:
#>  0
```
