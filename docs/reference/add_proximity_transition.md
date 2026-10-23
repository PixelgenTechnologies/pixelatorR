# Add a proximity transition column to differential proximity results

Classifies the reference/rest and target group medians into proximity
regimes around zero (`"anti"` for `< -tau`, `"random"` for
`[-tau, tau]`, and `"coloc"` for `> tau`) and combines them into a
`transition` column formatted as `"<ref_state> -> <tgt_state>"`.

## Usage

``` r
add_proximity_transition(
  data,
  tau = 0.15,
  tgt_col = "med_tgt",
  ref_col = NULL,
  transition_col = "transition",
  separator = " -> "
)
```

## Arguments

- data:

  A `data.frame` or `tbl_df` returned by
  [`DifferentialProximityAnalysis`](DifferentialProximityAnalysis.md) or
  [`FindAllProximityMarkers`](FindAllProximityMarkers.md).

- tau:

  Non-negative numeric threshold defining the neutral `"random"` zone
  `[-tau, tau]` around 0. Default is `0.15`.

- tgt_col:

  Name of the target median column. Default is `"med_tgt"`.

- ref_col:

  Name of the reference median column. When `NULL` (default),
  automatically uses `"med_ref"` or `"med_rest"` if present in `data`.

- transition_col:

  Name of the output transition column. Default is `"transition"`.

- separator:

  Separator placed between the reference and target states. Default is
  `" -> "`.

## Value

`data` with a character column named by `transition_col`. An existing
column of that name is replaced.

## Details

Works on output tables from both
[`DifferentialProximityAnalysis`](DifferentialProximityAnalysis.md)
(using `med_ref` and `med_tgt`) and
[`FindAllProximityMarkers`](FindAllProximityMarkers.md) (using
`med_rest` and `med_tgt`).

## Examples

``` r
library(dplyr)
df <- tibble(
  target = c("B", "CD4 T", "NK"),
  med_tgt = c(0.8, 0.9, -0.4),
  med_rest = c(0.0, 0.5, 0.0),
  diff_median = med_tgt - med_rest
)
add_proximity_transition(df, tau = 0.15)
#> # A tibble: 3 × 5
#>   target med_tgt med_rest diff_median transition     
#>   <chr>    <dbl>    <dbl>       <dbl> <chr>          
#> 1 B          0.8      0           0.8 random -> coloc
#> 2 CD4 T      0.9      0.5         0.4 coloc -> coloc 
#> 3 NK        -0.4      0          -0.4 random -> anti 
```
