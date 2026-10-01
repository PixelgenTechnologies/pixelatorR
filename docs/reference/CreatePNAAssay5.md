# Create a PNAAssay5 object

Create a [`PNAAssay5`](PNAAssay5-class.md) object from a count matrix
and a list of [`CellGraph`](CellGraph-class.md) objects. The expected
format of the input matrix is features x cells. Optionally, a `tbl_df`
with proximity scores and a `tbl_df` with information on source PXL file
paths can be provided.

## Usage

``` r
CreatePNAAssay5(
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
  [`CreateAssay5Object`](https://satijalab.github.io/seurat-object/reference/CreateAssay5Object.html)

## Value

A `PNAAssay5` object

## Examples

``` r
library(pixelatorR)
library(dplyr)

pxl_file <- minimal_pna_pxl_file()
counts <- ReadPNA_counts(pxl_file)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpmC3mql/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
pna_assay5 <- CreatePNAAssay5(
  counts = counts,
  cellgraphs = rep(list(NULL), ncol(counts)) %>%
    setNames(colnames(counts))
)
pna_assay5
#> PNAAssay (v5) data with 158 features for 5 cells
#> First 10 features:
#>  HLA-ABC, B2M, CD11b, CD11c, CD18, CD82, CD8, TCRab, HLA-DR, CD45 
#> Layers:
#>  counts 
#> Loaded CellGraph objects:
#>  0
```
