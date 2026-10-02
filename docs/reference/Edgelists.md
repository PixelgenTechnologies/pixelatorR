# Load edgelists

Get the edgelist(s) from a [`PNAAssay`](PNAAssay-class.md),
[`PNAAssay5`](PNAAssay5-class.md) or a `Seurat` object. The edgelist(s)
are stored on disk in the PXL files, so the method only works if the
paths are set correctly (see [`?FSMap`](FSMap.md)).

## Usage

``` r
Edgelists(object, ...)

# S3 method for class 'PNAAssay'
Edgelists(object, lazy = TRUE, union = TRUE, ...)

# S3 method for class 'PNAAssay5'
Edgelists(object, lazy = TRUE, union = TRUE, ...)

# S3 method for class 'Seurat'
Edgelists(
  object,
  assay = NULL,
  meta_data_columns = NULL,
  lazy = TRUE,
  union = TRUE,
  ...
)
```

## Arguments

- object:

  An object with polarization scores

- ...:

  Not implemented

- lazy:

  A logical indicating whether to lazy load the edgelist(s) from the PXL
  files

- union:

  A logical indicating whether to return the union of all edgelists from
  all PXL files (TRUE) or a list of edgelists (FALSE).

- assay:

  Name of a `CellGraphAssay`

- meta_data_columns:

  A character vector with meta.data column names. This option can be
  useful to join meta.data columns with the proximity score table.

## Value

`Edgelists`: Edgelist(s)

## See also

Other spatial metrics:
[`ColocalizationScores()`](ColocalizationScores.md),
[`PolarizationScores()`](PolarizationScores.md),
[`ProximityScores()`](ProximityScores.md)

## Examples

``` r
library(pixelatorR)

pxl_file <- minimal_pna_pxl_file()
seur_obj <- ReadPNA_Seurat(pxl_file)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpjKKHFf/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> ✔ Created a <Seurat> object with 5 cells and 158 targeted surface proteins
el <- Edgelists(seur_obj[["PNA"]], lazy = FALSE)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpjKKHFf/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
el
#> # A tibble: 528,594 × 7
#>    marker_1 marker_2    umi1    umi2 read_count uei_count component       
#>    <chr>    <chr>    <int64> <int64>    <int64>   <int64> <chr>           
#>  1 B2M      HLA-ABC    5 e15   2 e16          2         1 d4074c845bb62800
#>  2 B2M      HLA-ABC    3 e16   5.e16          2         1 d4074c845bb62800
#>  3 B2M      HLA-ABC    6 e16   6 e16          3         1 d4074c845bb62800
#>  4 B2M      HLA-ABC    5.e16   4 e16          1         1 d4074c845bb62800
#>  5 B2M      HLA-ABC    6 e16   6 e15          2         1 d4074c845bb62800
#>  6 B2M      HLA-ABC    5.e16   4 e16          1         1 d4074c845bb62800
#>  7 B2M      HLA-ABC    4 e16   6 e16          1         1 d4074c845bb62800
#>  8 B2M      HLA-ABC    6 e16   4 e15          1         1 d4074c845bb62800
#>  9 B2M      HLA-ABC    1.e16   1.e16          1         1 d4074c845bb62800
#> 10 B2M      HLA-ABC    1.e16   3 e16          1         1 d4074c845bb62800
#> # ℹ 528,584 more rows
```
