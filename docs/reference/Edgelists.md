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
#> ℹ /tmp/Rtmp3BjCtL/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> ✔ Created a <Seurat> object with 5 cells and 158 targeted surface proteins
el <- Edgelists(seur_obj[["PNA"]], lazy = FALSE)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/Rtmp3BjCtL/duckdb
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
#>  1 CD102    CD45RA     1.e16   4 e16          1         1 2708240b908e2eba
#>  2 CD102    CD45RA     1.e16   7 e16          2         1 2708240b908e2eba
#>  3 CD102    CD45RA     5 e15   1 e16          2         2 2708240b908e2eba
#>  4 CD102    CD45RA     1.e16   1 e16          2         1 2708240b908e2eba
#>  5 CD102    CD45RA     1.e16   5.e16          2         1 2708240b908e2eba
#>  6 CD102    CD45RA     7 e15   6 e16          1         1 2708240b908e2eba
#>  7 CD102    CD45RA     2 e15   6 e16          1         1 2708240b908e2eba
#>  8 CD102    CD45RA     3 e16   3 e16          3         1 2708240b908e2eba
#>  9 CD102    CD45RA     5 e15   6 e16          1         1 2708240b908e2eba
#> 10 CD102    CD45RA     4 e16   5.e16          1         1 2708240b908e2eba
#> # ℹ 528,584 more rows
```
