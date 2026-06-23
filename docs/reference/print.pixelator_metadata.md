# Print method for pixelator_metadata

Print method for pixelator_metadata

## Usage

``` r
# S3 method for class 'pixelator_metadata'
print(x, detailed = TRUE, ...)
```

## Arguments

- x:

  An object of class `pixelator_metadata`

- detailed:

  Logical. If `TRUE` print all metadata, if `FALSE` print only the
  pixelator version and sample ID.

- ...:

  Additional arguments passed to
  [`print`](https://satijalab.github.io/seurat-object/reference/print.DimReduc.html)

## Examples

``` r
library(pixelatorR)
library(dplyr)

# Load example data
pxl_file <- minimal_mpx_pxl_file()
meta_data <- ReadMPX_metadata(pxl_file)

# Check pixelator version and sample ID
meta_data
#> Sample 1 name:       Mock_data
#> Pixelator version:   0.12.0
#> 
#> 
#> ── Analysis parameters ──
#> 
#> # A tibble: 9 × 2
#>   parameter                         value
#>   <chr>                             <chr>
#> 1 compute_polarization              TRUE 
#> 2 compute_colocalization            TRUE 
#> 3 use_full_bipartite                FALSE
#> 4 polarization_normalization        clr  
#> 5 polarization_binarization         FALSE
#> 6 colocalization_transformation     log1p
#> 7 colocalization_neighbourhood_size 1    
#> 8 colocalization_n_permutations     50   
#> 9 colocalization_min_region_count   5    

# Multiple files with less detail
pxl_files <- c(pxl_file, pxl_file)
meta_data_merged <- lapply(pxl_files, ReadMPX_metadata) %>%
  bind_rows()
meta_data_merged %>% print(detailed = FALSE)
#> Sample 1 name:       Mock_data
#> Pixelator version:   0.12.0
#> 
#> ────────────────────────────────────────────────────────────────────────────────
#> Sample 2 name:       Mock_data
#> Pixelator version:   0.12.0
#> 
```
