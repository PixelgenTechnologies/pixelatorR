# Calculate the sizes of the LCCs from an edgelist

This function calculates the sizes of the largest connected components
(LCCs) from a set of parquet files containing downsampled edgelists. The
LCCs are computed for each fraction specified in the input tibble `df`,
linking each fraction to its corresponding parquet file. Moreover, the
LCCs are computed per cell (component).

## Usage

``` r
lcc_sizes(df, mc_cores = 1)
```

## Arguments

- df:

  A tibble with columns `pq_files` (paths to parquet files) and `fracs`
  generated with [`downsample_to_parquet`](downsample_to_parquet.md).

- mc_cores:

  Number of cores to use for parallel processing.

## Value

A tibble with the sizes of the largest connected components (LCC) for
each fraction in `df$fracs`. The results are grouped by cell
(component).
