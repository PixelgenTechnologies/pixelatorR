# Summarize proximity scores

Computes the median or mean proximity score for each protein pair across
all group combinations defined by `group_vars`. This is typically useful
when profiling a population of interest, where the goal is to compute
representative proximity scores for each protein pair within that
population.

The proximity table typically contains a small fraction of the possible
protein pairs for each cell. These missing pairs typically have UMI
counts below the detection threshold and their proximity scores can
therefore be imputed as 0 (0 observed join counts and 0 expected join
counts =\> no deviation from the expected value).
`SummarizeProximityScores` pads the proximity score vectors with 0s for
the missing pairs to ensure that the summary statistics are computed for
the entire population. This behavior can be turned off by setting
`include_missing_obs = FALSE`.

## Usage

``` r
SummarizeProximityScores(object, ...)

# S3 method for class 'tbl_lazy'
SummarizeProximityScores(
  object,
  proximity_metric = c("log2_ratio", "join_count_z"),
  group_vars = NULL,
  include_missing_obs = TRUE,
  summary_stat = c("mean", "median"),
  detailed = FALSE,
  ...
)

# S3 method for class 'data.frame'
SummarizeProximityScores(
  object,
  proximity_metric = "log2_ratio",
  group_vars = NULL,
  include_missing_obs = TRUE,
  summary_stat = c("mean", "median"),
  detailed = FALSE,
  ...
)
```

## Arguments

- object:

  A `tbl_df` or `tbl_lazy` object with proximity scores.

- ...:

  Additional arguments. Currently not used.

- proximity_metric:

  The proximity metric to use. One of "log2_ratio" or "join_count_z".

- group_vars:

  A character vector with the names of the variables to use for grouping
  This is typically used is you want to summarize the proximity scores
  for multiple cell types and/or conditions.

- include_missing_obs:

  Logical indicating whether to include missing observations as 0 when
  computing the summary statistics.

- summary_stat:

  One of "mean" or "median"

- detailed:

  Logical indicating whether to return lists which can be used to
  compute custom summary statistics. See examples for details

## Value

A `tbl_df` with summary statistics

## Examples

``` r
library(pixelatorR)
library(dplyr)

pxl_file <- minimal_pna_pxl_file()
se <- ReadPNA_Seurat(pxl_file)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpjKKHFf/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> ✔ Created a <Seurat> object with 5 cells and 158 targeted surface proteins
proximity_table <- ProximityScores(se)

# Default method uses mean
SummarizeProximityScores(proximity_table) %>% head()
#> # A tibble: 6 × 7
#>   marker_1 marker_2 n_cells_detected n_cells n_cells_missing pct_detected
#>   <chr>    <chr>               <int>   <int>           <int>        <dbl>
#> 1 CD56     Siglec-9                4       5               1          0.8
#> 2 CD56     CD79a                   4       5               1          0.8
#> 3 CD56     HLA-ABC                 4       5               1          0.8
#> 4 CD56     GPR56                   3       5               2          0.6
#> 5 CD56     CD82                    4       5               1          0.8
#> 6 CD56     IgD                     4       5               1          0.8
#> # ℹ 1 more variable: mean_log2_ratio <dbl>

# Switch to median
SummarizeProximityScores(proximity_table, summary_stat = "median") %>% head()
#> # A tibble: 6 × 7
#>   marker_1 marker_2 n_cells_detected n_cells n_cells_missing pct_detected
#>   <chr>    <chr>               <int>   <int>           <int>        <dbl>
#> 1 CD56     Siglec-9                4       5               1          0.8
#> 2 CD56     CD79a                   4       5               1          0.8
#> 3 CD56     HLA-ABC                 4       5               1          0.8
#> 4 CD56     GPR56                   3       5               2          0.6
#> 5 CD56     CD82                    4       5               1          0.8
#> 6 CD56     IgD                     4       5               1          0.8
#> # ℹ 1 more variable: median_log2_ratio <dbl>

# Ignore missing values
SummarizeProximityScores(proximity_table, include_missing_obs = FALSE) %>% head()
#> # A tibble: 6 × 7
#>   marker_1 marker_2 n_cells_detected n_cells n_cells_missing pct_detected
#>   <chr>    <chr>               <int>   <int>           <int>        <dbl>
#> 1 CD56     Siglec-9                4       5               1          0.8
#> 2 CD56     CD79a                   4       5               1          0.8
#> 3 CD56     HLA-ABC                 4       5               1          0.8
#> 4 CD56     GPR56                   3       5               2          0.6
#> 5 CD56     CD82                    4       5               1          0.8
#> 6 CD56     IgD                     4       5               1          0.8
#> # ℹ 1 more variable: mean_log2_ratio <dbl>

# Return lists which can be used to compute custom summary statistics
SummarizeProximityScores(proximity_table, detailed = TRUE) %>%
  # It's important to do rowwise computations
  rowwise() %>%
  mutate(
    sd = sd(unlist(log2_ratio_list)),
    iqr = IQR(unlist(log2_ratio_list)),
    mad = mad(unlist(log2_ratio_list)),
    q90 = quantile(unlist(log2_ratio_list), 0.9)
  ) %>%
  select(marker_1, marker_2, sd, iqr, mad, q90) %>%
  ungroup()
#> # A tibble: 12,561 × 6
#>    marker_1 marker_2    sd   iqr   mad   q90
#>    <chr>    <chr>    <dbl> <dbl> <dbl> <dbl>
#>  1 CD56     CD71     0     0         0 0    
#>  2 CD56     IgM      0     0         0 0    
#>  3 CD56     TCRva7.2 0     0         0 0    
#>  4 CD56     KLRG1    0     0         0 0    
#>  5 CD56     CD93     0     0         0 0    
#>  6 CD56     HLA-DR   0.128 0         0 0    
#>  7 CD56     CD69     0     0         0 0    
#>  8 CD56     CD59     0.826 0.258     0 0    
#>  9 CD56     CD62P    0     0         0 0    
#> 10 CD56     CX3CR1   0.709 0         0 0.951
#> # ℹ 12,551 more rows
```
