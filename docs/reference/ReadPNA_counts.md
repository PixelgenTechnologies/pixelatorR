# Read a count matrix from a PXL file with PNA data

Read a count matrix from a PXL file with PNA data

## Usage

``` r
ReadPNA_counts(pxl_file)
```

## Arguments

- pxl_file:

  Path to a PXL file

## Value

A count matrix or a list if `return_list = TRUE`

## See also

Other PXL-data-loaders: [`ReadPNA_Seurat()`](ReadPNA_Seurat.md),
[`ReadPNA_metadata()`](ReadPNA_metadata.md)

## Examples

``` r
library(pixelatorR)

# Load example counts
pxl_file <- minimal_pna_pxl_file()
counts <- ReadPNA_counts(pxl_file)
counts[1:5, 1:5]
#> 5 x 5 sparse Matrix of class "dgCMatrix"
#>         0a45497c6bfbfb22 2708240b908e2eba c3c393e9a17c1981 d4074c845bb62800
#> HLA-ABC              865             2077             2480             2994
#> B2M                 1182             3448             5307             9753
#> CD11b                929               12               17                6
#> CD11c                133             1354               26                8
#> CD18                 716              780              326             1127
#>         efe0ed189cb499fc
#> HLA-ABC             2212
#> B2M                 6082
#> CD11b                  .
#> CD11c                  4
#> CD18                 523
```
