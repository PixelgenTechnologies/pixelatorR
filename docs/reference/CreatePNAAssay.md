# Create a PNAAssay object

Create a [`PNAAssay`](PNAAssay-class.md) object from a count matrix and
a list of [`CellGraph`](CellGraph-class.md) objects. The expected format
of the input matrix is features x cells. Optionally, a `tbl_df` with
proximity scores and a `tbl_df` with information on source PXL file
paths can be provided.

## Usage

``` r
CreatePNAAssay(
  counts,
  cellgraphs,
  proximity = NULL,
  fs_map = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- counts:

  Unnormalized data (raw counts)

- cellgraphs:

  A named list with [`CellGraph`](CellGraph-class.md) objects

- proximity:

  A `tbl_df` with proximity scores

- fs_map:

  A `tbl_df` with information on source PXL file paths, sample IDs, and
  component IDs

- verbose:

  Print messages

- ...:

  Additional arguments passed to
  [`CreateAssayObject`](https://satijalab.github.io/seurat-object/reference/CreateAssayObject.html)

## Value

A `PNAAssay` object

## Examples

``` r
library(pixelatorR)
library(dplyr)

pxl_file <- minimal_pna_pxl_file()
counts <- ReadPNA_counts(pxl_file)
pna_assay <- CreatePNAAssay(
  counts = counts,
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
