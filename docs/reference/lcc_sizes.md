# Compute LCC sizes from downsampled edgelists

Calculates the sizes of the largest connected components (LCCs) from a
set of parquet files containing downsampled edgelists.

## Usage

``` r
lcc_sizes(df, mc_cores = 1)
```

## Arguments

- df:

  A tibble with columns:

  pq_files

  :   Paths to parquet files containing downsampled edgelists.

  fracs

  :   The corresponding downsampling fractions.

  Typically generated with
  [`downsample_to_parquet`](downsample_to_parquet.md).

- mc_cores:

  Number of cores to use for parallel processing. Default is 1
  (sequential processing). Values \> 1 use
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html).

## Value

A tibble with columns:

- component:

  The component (cell) identifier.

- frac:

  The downsampling fraction.

- n_nodes:

  The number of nodes in the largest connected component.

## Details

### Algorithm

For each parquet file (representing a downsampling fraction), this
function:

1.  Loads the edgelist into a DuckDB connection

2.  Constructs a property graph using the `duckpgq` extension

3.  Computes weakly connected components

4.  Identifies the largest component for each cell

### Performance

The `duckpgq` extension enables efficient graph operations without
loading the full edgelist into memory. Parallel processing across
fractions is supported via `mc_cores`.

## See also

[`downsample_to_parquet`](downsample_to_parquet.md) for generating the
input parquet files. [`lcc_curve`](lcc_curve.md) for a high-level
wrapper.
