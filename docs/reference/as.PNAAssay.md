# Convert objects to a [`PNAAssay`](PNAAssay-class.md)

Convert objects to a [`PNAAssay`](PNAAssay-class.md)

## Usage

``` r
as.PNAAssay(x, ...)

# S3 method for class 'Assay'
as.PNAAssay(x, cellgraphs = NULL, proximity = NULL, fs_map = NULL, ...)
```

## Arguments

- x:

  An object to convert to class [`PNAAssay`](PNAAssay-class.md)

- ...:

  Arguments passed to other methods

- cellgraphs:

  A list of [`CellGraph`](CellGraph-class.md) objects

- proximity:

  A `tbl_df` with proximity scores

- fs_map:

  A `tbl_df` with information on source PXL file paths, sample IDs, and
  component IDs

## Value

A `PNAAssay` object

## Examples

``` r
library(pixelatorR)
library(dplyr)
library(SeuratObject)

pxl_file <- minimal_pna_pxl_file()
counts <- ReadPNA_counts(pxl_file)
assay <- CreateAssayObject(
  counts = counts
)
pna_assay <- as.PNAAssay(
  assay,
  cellgraphs = rep(list(NULL), ncol(counts)) %>%
    setNames(colnames(counts))
)
pna_assay
#> PNAAssay data with 158 features for 5 cells
#> First 10 features:
#>  HLA-ABC, B2M, CD11b, CD11c, CD18, CD82, CD8, TCRab, HLA-DR, CD45 
#> Loaded CellGraph objects:
#>  0
```
