# ProximityScores

Get and set proximity scores for a [`PNAAssay`](PNAAssay-class.md),
[`PNAAssay5`](PNAAssay5-class.md) or a `Seurat` object

## Usage

``` r
ProximityScores(object, ...)

ProximityScores(object, ...) <- value

# S3 method for class 'PNAAssay'
ProximityScores(
  object,
  add_marker_counts = FALSE,
  add_marker_proportions = FALSE,
  lazy = FALSE,
  calc_log2ratio = TRUE,
  ...
)

# S3 method for class 'PNAAssay5'
ProximityScores(
  object,
  add_marker_counts = FALSE,
  add_marker_proportions = FALSE,
  lazy = FALSE,
  calc_log2ratio = TRUE,
  ...
)

# S3 method for class 'PNAAssay'
ProximityScores(object, ...) <- value

# S3 method for class 'PNAAssay5'
ProximityScores(object, ...) <- value

# S3 method for class 'Seurat'
ProximityScores(
  object,
  assay = NULL,
  meta_data_columns = NULL,
  add_marker_counts = FALSE,
  add_marker_proportions = FALSE,
  lazy = FALSE,
  calc_log2ratio = TRUE,
  ...
)

# S3 method for class 'Seurat'
ProximityScores(object, assay = NULL, ...) <- value
```

## Arguments

- object:

  An object with polarization scores

- ...:

  Not implemented

- value:

  A `tbl_df` or `tbl_lazy` with proximity scores

- add_marker_counts:

  A logical indicating whether to add marker counts to the output
  ("count_1" and "count_2")

- add_marker_proportions:

  A logical indicating whether to add marker count proportions to the
  output ("p1" and "p2")

- lazy:

  A logical indicating whether to lazy load the proximity scores from
  the PXL files

- calc_log2ratio:

  A logical indicating whether to calculate the log2 ratio proximity
  score

- assay:

  Name of a `PNAAssay`

- meta_data_columns:

  A character vector with meta.data column names. This option can be
  useful to join meta.data columns with the proximity score table.

## Value

`ProximityScores`: Proximity scores

`ProximityScores<-`: An object with updated proximity scores

## See also

Other spatial metrics:
[`ColocalizationScores()`](ColocalizationScores.md),
[`Edgelists()`](Edgelists.md),
[`PolarizationScores()`](PolarizationScores.md)

## Examples

``` r
library(pixelatorR)
library(dplyr)

pxl_file <- minimal_pna_pxl_file()
seur_obj <- ReadPNA_Seurat(pxl_file)
ProximityScores(seur_obj[["PNA"]])
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

# Set proximity scores
ProximityScores(seur_obj[["PNA"]]) <-
  ProximityScores(seur_obj[["PNA"]]) %>%
  mutate(ratio = join_count / join_count_expected_mean)

library(pixelatorR)

# Create example Seurat object
pxl_file <- minimal_pna_pxl_file()
seur_obj <- ReadPNA_Seurat(pxl_file)

# Get proximity scores
proximity <- ProximityScores(seur_obj)
proximity
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

# Get proximity scores with additional meta data
proximity <-
  ProximityScores(seur_obj, meta_data_columns = c("n_umi", "n_edges"))
proximity
#> # A tibble: 58,696 × 11
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
#> # ℹ 6 more variables: join_count_z <dbl>, join_count_p <dbl>, component <chr>,
#> #   log2_ratio <dbl>, n_umi <int>, n_edges <int>

# Get proximity scores with marker_1 and marker_2 counts
proximity <-
  ProximityScores(seur_obj, add_marker_counts = TRUE)
proximity
#> # A tibble: 58,696 × 11
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
#> # ℹ 6 more variables: join_count_z <dbl>, join_count_p <dbl>, component <chr>,
#> #   log2_ratio <dbl>, count_1 <int>, count_2 <int>

# Update proximity scores in Seurat object
ProximityScores(seur_obj) <- ProximityScores(seur_obj) %>%
  mutate(ratio = join_count / join_count_expected_mean)
```
