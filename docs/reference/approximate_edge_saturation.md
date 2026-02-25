# Compute approximate edge saturation

This function computes an approximate edge saturation for each component
in the edgelist. It estimates the theoretical maximum number of edges
using the Chao1 estimator and calculates the edge saturation as the
ratio of the actual number of edges to the theoretical maximum.

## Usage

``` r
approximate_edge_saturation(db, components = NULL, table_name = NULL)
```

## Arguments

- db:

  A `PixelDB` object.

- components:

  Optional character vector of component names to filter the edgelist.

- table_name:

  Optional name of the remote table.

## Value

A `lazy_df` with the edge saturation.
