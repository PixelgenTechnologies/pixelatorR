# Read metadata from a PXL file

Read metadata from a PXL file

## Usage

``` r
ReadMPX_metadata(filename)
```

## Arguments

- filename:

  Path to a PXL file

## Examples

``` r
library(pixelatorR)

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

# Check parameter settings for the pixelator run
meta_data$analysis$params
#>              compute_polarization            compute_colocalization 
#>                            "TRUE"                            "TRUE" 
#>                use_full_bipartite        polarization_normalization 
#>                           "FALSE"                             "clr" 
#>         polarization_binarization     colocalization_transformation 
#>                           "FALSE"                           "log1p" 
#> colocalization_neighbourhood_size     colocalization_n_permutations 
#>                               "1"                              "50" 
#>   colocalization_min_region_count 
#>                               "5" 
```
