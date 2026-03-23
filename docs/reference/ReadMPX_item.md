# Read an MPX data item

`ReadMPX_item` reads any number of items from a PXL file. If multiple
items are specified, the output is a list of `tbl_df` objects. Otherwise
the output is a single `tbl_df`. `ReadMPX_polarization`,
`ReadMPX_colocalization` and `ReadMPX_edgelist` are wrappers for
`ReadMPX_item` to read specific items.

## Usage

``` r
ReadMPX_item(
  filename,
  items = c("colocalization", "polarization", "edgelist"),
  verbose = TRUE
)

ReadMPX_polarization(filename, verbose = TRUE)

ReadMPX_colocalization(filename, verbose = TRUE)

ReadMPX_edgelist(filename, verbose = TRUE)
```

## Arguments

- filename:

  Path to a PXL file

- items:

  One or several of "colocalization", "polarization", "edgelist"

- verbose:

  Print messages

## Value

Either an object of class `tbl_df` or a list of `tbl_df` objects if
multiple `items` are selected

## See also

Other data-loaders: [`ReadMPX_Seurat()`](ReadMPX_Seurat.md),
[`ReadPNA_proximity()`](ReadPNA_proximity.md)

## Examples

``` r
library(pixelatorR)

# Load example data
pxl_file <- minimal_mpx_pxl_file()
polarization <- ReadMPX_item(pxl_file, items = "polarization")
#> ! Failed to remove temporary dir C:/Users/max/AppData/Local/Temp/Rtmpampkmn/dir5bf42d6069d3
polarization
#> # A tibble: 380 × 6
#>    morans_i morans_p_value morans_p_adjusted morans_z marker component    
#>       <dbl>          <dbl>             <dbl>    <dbl> <chr>  <chr>        
#>  1 -0.00165          0.894             0.999   -0.134 ACTB   RCVCMP0000217
#>  2 -0.00734          0.499             0.999   -0.676 B2M    RCVCMP0000217
#>  3 -0.00850          0.409             0.999   -0.826 CD102  RCVCMP0000217
#>  4  0.00102          0.856             0.999    0.182 CD11a  RCVCMP0000217
#>  5 -0.00381          0.656             0.999   -0.445 CD11b  RCVCMP0000217
#>  6 -0.00493          0.631             0.999   -0.480 CD11c  RCVCMP0000217
#>  7 -0.00497          0.596             0.999   -0.530 CD127  RCVCMP0000217
#>  8 -0.00150          0.910             0.999   -0.113 CD137  RCVCMP0000217
#>  9 -0.00224          0.850             0.999   -0.190 CD14   RCVCMP0000217
#> 10 -0.00246          0.837             0.999   -0.205 CD150  RCVCMP0000217
#> # ℹ 370 more rows

# Alternative 2
polarization <- ReadMPX_polarization(pxl_file)
#> ! Failed to remove temporary dir C:/Users/max/AppData/Local/Temp/Rtmpampkmn/dir5bf46b0270d2
```
