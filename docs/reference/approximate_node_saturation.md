# Compute approximate node saturation

This function computes an approximate node saturation for each component
in the edgelist. It estimates the theoretical maximum number of nodes
using the Chao1 estimator and calculates the nodes saturation as the
ratio of the actual number of nodes to the theoretical maximum.

## Usage

``` r
approximate_node_saturation(db, components = NULL, table_name = NULL)
```

## Arguments

- db:

  A `PixelDB` object.

- components:

  Optional vector of component names to filter the edgelist by.

- table_name:

  Optional name of the remote table.

## Value

A `lazy_df` with the node saturation.
