# Load the Proximity scores table from a PNA PXL file

Load the Proximity scores table from a PNA PXL file

## Usage

``` r
ReadPNA_proximity(
  pxl_file,
  calc_log2_ratio = TRUE,
  lazy = FALSE,
  verbose = TRUE
)
```

## Arguments

- pxl_file:

  Path to a PXL file containing PNA data

- calc_log2_ratio:

  A logical specifying whether to calculate and add a log2ratio column
  to the output table. Default is `TRUE`

- lazy:

  A logical specifying whether to load the data lazily. If `TRUE`, a
  `tbl_lazy` object is returned.

- verbose:

  Print messages

## Value

A `tbl_df` or a `tbl_lazy` with PNA Proximity scores

## See also

Other data-loaders: [`ReadMPX_Seurat()`](ReadMPX_Seurat.md),
[`ReadMPX_item()`](ReadMPX_item.md)

## Examples

``` r
library(pixelatorR)

pxl_file <- minimal_pna_pxl_file()
proximity_tbl <- ReadPNA_proximity(pxl_file)
proximity_tbl
#> # A tibble: 58,696 × 9
#>    marker_1 marker_2 join_count join_count_expected_mean join_count_expected_sd
#>    <chr>    <chr>         <dbl>                    <dbl>                  <dbl>
#>  1 CD56     CD56              0                     0                     0    
#>  2 CD56     mIgG2b            0                     0.03                  0.171
#>  3 CD56     CD71              0                     0.01                  0.1  
#>  4 CD56     CD6               0                     2.07                  1.58 
#>  5 CD56     Siglec-9          0                     0.03                  0.171
#>  6 CD56     CD79a             0                     0                     0    
#>  7 CD56     NKp80             0                     0                     0    
#>  8 CD56     CD85j             0                     0.17                  0.378
#>  9 CD56     IgM               0                     0.08                  0.307
#> 10 CD56     TCRva7.2          0                     0.01                  0.1  
#> # ℹ 58,686 more rows
#> # ℹ 4 more variables: join_count_z <dbl>, join_count_p <dbl>, component <chr>,
#> #   log2_ratio <dbl>
```
