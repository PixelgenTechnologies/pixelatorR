# Read a count matrix from a pxl file

Read a count matrix from a pxl file

## Usage

``` r
ReadMPX_counts(filename, return_list = FALSE, verbose = TRUE)
```

## Arguments

- filename:

  Path to a PXL file

- return_list:

  If TRUE, returns a list with the expression matrix and the output from
  `h5read`

- verbose:

  Print messages

## Value

A count matrix or a list if `return_list = TRUE`

## Examples

``` r
library(pixelatorR)

# Load example data
pxl_file <- minimal_mpx_pxl_file()
counts <- ReadMPX_counts(pxl_file)
counts[1:5, 1:5]
#>       RCVCMP0000217 RCVCMP0000118 RCVCMP0000487 RCVCMP0000655 RCVCMP0000263
#> CD274            18             6             9            23            22
#> CD44            117            49           100           338           417
#> CD25             14             5            13             9            34
#> CD279            39             1             7             8             2
#> CD41              2            18            18             9             2
```
