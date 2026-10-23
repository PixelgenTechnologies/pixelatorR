# Find all proximity markers

Identifies proximity marker pairs that differ between each group and the
cells in every other group. This is the proximity-score equivalent of
`FindAllMarkers` in Seurat.

## Usage

``` r
FindAllProximityMarkers(object, ...)

# S3 method for class 'Seurat'
FindAllProximityMarkers(
  object,
  group_by,
  idents = NULL,
  assay = NULL,
  lazy = FALSE,
  min_exp_join_count = 0,
  min_cells_per_group = 10,
  diff_threshold = 0.01,
  min_pct = 0,
  min_diff_pct = -Inf,
  only_pos = FALSE,
  proximity_metric = "log2_ratio",
  metric_type = c("all", "self", "co"),
  p_adjust_method = c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr"),
  verbose = TRUE,
  ...
)

# S3 method for class 'data.frame'
FindAllProximityMarkers(
  object,
  group_by,
  idents = NULL,
  min_cells_per_group = 10,
  diff_threshold = 0.01,
  min_pct = 0,
  min_diff_pct = -Inf,
  only_pos = FALSE,
  proximity_metric = "log2_ratio",
  metric_type = c("all", "self", "co"),
  p_adjust_method = c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr"),
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  An object containing PNA proximity scores.

- ...:

  Not used. `contrast_column`, `reference`, `targets`, `group_vars`, and
  `group_data` are set by this function and cannot be supplied. `method`
  and `backend` are not used.

- group_by:

  Name of a single column of group labels. The column must be a
  character vector or a factor with at least two levels.

- idents:

  Optional levels of `group_by` to test. Numeric and factor values are
  coerced to character. When `NULL` (default), every level is tested
  against the others.

- assay:

  Name of the assay to use.

- lazy:

  If `TRUE`, proximity scores are fetched with
  [`ProximityScores`](ProximityScores.md)`(object, lazy = TRUE)` and
  filtered with the `duckdb` backend before testing. Set this when the
  `Seurat` object was created with `load_proximity_scores = FALSE` and
  the proximity scores were not stored in the object.

- min_exp_join_count:

  Minimum expected join count for a marker pair to be included. Pairs
  below the threshold are treated as missing and set to 0.

- min_cells_per_group:

  Minimum number of cells required in both the level being tested and
  the pooled remainder. Comparisons that do not meet this threshold are
  skipped. The default is 10.

- diff_threshold:

  Minimum absolute difference in the proximity metric required to test a
  marker pair. The default is 0.01. When `only_pos` is `TRUE`, the
  difference must also be positive.

- min_pct:

  Minimum fraction of cells in either group with a non-zero score. The
  default is 0.

- min_diff_pct:

  Minimum difference in the fraction of cells with a non-zero score. The
  default is `-Inf`.

- only_pos:

  If `TRUE`, test only marker pairs whose median is higher in the level
  being tested than in the pooled remainder. The default is `FALSE`.

- proximity_metric:

  The proximity metric to use. Any numeric column in the proximity score
  table can be selected. The default is `"log2_ratio"`.

- metric_type:

  One of `"all"`, `"self"`, or `"co"`. `"self"` keeps pairs where both
  markers are the same protein, and `"co"` keeps pairs of different
  proteins.

- p_adjust_method:

  One of `"bonferroni"`, `"holm"`, `"hochberg"`, `"hommel"`, `"BH"`,
  `"BY"`, or `"fdr"` (see
  [`?p.adjust`](https://rdrr.io/r/stats/p.adjust.html)). P-values are
  adjusted within each one-versus-rest comparison.

- verbose:

  Print messages.

## Value

A `tbl_df` of test results stacked across comparisons.

## Details

`group_by` must be a single column, for example `"seurat_clusters"` or
`"cell_type"`. For each level of that column, the level is compared to
the remaining levels pooled into one reference group. Each comparison is
delegated to
[`DifferentialProximityAnalysis`](DifferentialProximityAnalysis.md).
With clusters `0`, `1`, and `2`, the comparisons are `0` versus the pool
of `1` and `2`, `1` versus the pool of `0` and `2`, and `2` versus the
pool of `0` and `1`.

The level being tested is stored in `target`. The pooled remainder is
stored in `reference` and is labeled `"rest"`, unless that label is
already a level of `group_by`. A positive `diff_median` means the
proximity metric is higher in `target` than in the pooled remainder. Set
`only_pos = TRUE` to test only those pairs. Factor levels set the
comparison order. Character labels are tested in order of first
appearance. P-values are adjusted separately within each comparison,
across the pairs that were tested.

Proximity scores are pivoted once with
[`ProximityScoresToAssay`](ProximityScoresToAssay.md). Each comparison
then calls the matrix method of
[`DifferentialProximityAnalysis`](DifferentialProximityAnalysis.md) on
that same matrix. Missing scores are treated as 0.

A `Seurat` object created with `load_proximity_scores = FALSE` does not
store proximity scores. Pass `lazy = TRUE` to fetch them with
[`ProximityScores`](ProximityScores.md)`(object, lazy = TRUE)`.

## See also

Other DA-methods:
[`DifferentialProximityAnalysis()`](DifferentialProximityAnalysis.md),
[`RunDAA()`](RunDAA.md), [`RunDCA()`](RunDCA.md),
[`RunDPA()`](RunDPA.md)

## Examples

``` r
library(dplyr)
set.seed(1)
example_data <- tidyr::expand_grid(
  marker_1 = c("CD19", "CD3E"),
  marker_2 = c("CD19", "CD3E"),
  component = paste0("cell", seq_len(12))
) %>%
  mutate(
    join_count_z = rnorm(n()),
    seurat_clusters = dplyr::case_when(
      component %in% paste0("cell", 1:4) ~ "0",
      component %in% paste0("cell", 5:8) ~ "1",
      TRUE ~ "2"
    )
  )

# Compare each cluster with the cells in every other cluster
FindAllProximityMarkers(
  example_data,
  group_by = "seurat_clusters",
  proximity_metric = "join_count_z",
  min_cells_per_group = 4,
  verbose = FALSE
)
#> # A tibble: 12 × 10
#>    target pct_tgt pct_ref med_tgt med_rest diff_median     p p_adj marker_1
#>    <chr>    <dbl>   <dbl>   <dbl>    <dbl>       <dbl> <dbl> <dbl> <chr>   
#>  1 0            1       1  -0.221   0.439      -0.660  0.445     1 CD19    
#>  2 0            1       1  -0.333   0.688      -1.02   0.350     1 CD19    
#>  3 0            1       1  -0.106  -0.0783     -0.0277 0.799     1 CD3E    
#>  4 0            1       1   0.352   0.100       0.252  0.445     1 CD3E    
#>  5 1            1       1   0.408   0.287       0.122  1         1 CD19    
#>  6 1            1       1   0.708   0.0148      0.693  0.350     1 CD19    
#>  7 1            1       1   0.158  -0.106       0.264  0.445     1 CD3E    
#>  8 1            1       1   0.196   0.153       0.0434 1         1 CD3E    
#>  9 2            1       1   0.483   0.257       0.226  0.445     1 CD19    
#> 10 2            1       1   0.428   0.289       0.139  1         1 CD19    
#> 11 2            1       1  -0.234  -0.0795     -0.155  0.671     1 CD3E    
#> 12 2            1       1  -0.162   0.249      -0.411  0.445     1 CD3E    
#> # ℹ 1 more variable: marker_2 <chr>
```
