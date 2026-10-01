# Read metadata from a PNA PXL file

Read metadata from a PNA PXL file

## Usage

``` r
ReadPNA_metadata(pxl_file)
```

## Arguments

- pxl_file:

  Path to a PXL file

## Value

A `tbl_df` with sample meta data

## See also

Other PXL-data-loaders: [`ReadPNA_Seurat()`](ReadPNA_Seurat.md),
[`ReadPNA_counts()`](ReadPNA_counts.md)

## Examples

``` r
library(pixelatorR)

# Create example Seurat object
pxl_file <- minimal_pna_pxl_file()
ReadPNA_metadata(pxl_file)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpjKKHFf/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> # A tibble: 1 × 6
#>   sample_name        version technology   panel_name panel_version post_analysis
#>   <chr>              <chr>   <chr>        <chr>      <chr>         <named list> 
#> 1 PNA055_Sample07_S7 0.1.0   single-cell… pna-rnd-1… 0.1.0         <named list> 
```
