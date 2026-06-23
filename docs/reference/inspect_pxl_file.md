# Inspect a PXL file with MPX data

This function inspects a PXL file and returns a tibble with information
about the files contained in the PXL file.

## Usage

``` r
inspect_pxl_file(pxl_file)
```

## Arguments

- pxl_file:

  Path to a PXL file with MPX data

## Value

A tibble with information about the files contained in the PXL file

## Examples

``` r
pxl_file <- minimal_mpx_pxl_file()
inspect_pxl_file(pxl_file)
#> # A tibble: 5 × 3
#>   file_type                  n file     
#>   <chr>                  <int> <list>   
#> 1 adata.h5ad                 1 <chr [1]>
#> 2 edgelist.parquet           1 <chr [1]>
#> 3 metadata.json              1 <chr [1]>
#> 4 polarization.parquet       1 <chr [1]>
#> 5 colocalization.parquet     1 <chr [1]>
```
