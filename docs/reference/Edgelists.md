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
#> ✔ Created a <Seurat> object with 5 cells and 158 targeted surface proteins
el <- Edgelists(seur_obj[["PNA"]], lazy = FALSE)
el
#> # A tibble: 528,594 × 7
#>    marker_1 marker_2    umi1    umi2 read_count uei_count component       
#>    <chr>    <chr>    <int64> <int64>    <int64>   <int64> <chr>           
#>  1 CD6      B2M        5.e16   8 e15          1         1 c3c393e9a17c1981
#>  2 CD6      B2M        2.e16   4 e16          3         3 c3c393e9a17c1981
#>  3 CD6      B2M        4 e16   6 e16          1         1 c3c393e9a17c1981
#>  4 CD6      B2M        6 e16   6 e16          1         1 c3c393e9a17c1981
#>  5 CD6      B2M        2.e16   6 e16          1         1 c3c393e9a17c1981
#>  6 CD6      B2M        1.e16   6 e16          1         1 c3c393e9a17c1981
#>  7 CD6      B2M        3 e16   5.e16          2         2 c3c393e9a17c1981
#>  8 CD6      B2M        6 e16   4 e16          1         1 c3c393e9a17c1981
#>  9 CD6      B2M        6 e16   5.e16          5         3 c3c393e9a17c1981
#> 10 CD6      B2M        4 e16   5.e16          1         1 c3c393e9a17c1981
#> # ℹ 528,584 more rows
```
