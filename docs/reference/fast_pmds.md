# Fast pMDS implementation using RSpectra

This function computes a pMDS layout using the RSpectra package for
efficient eigen decomposition.

## Usage

``` r
fast_pmds(g, pivots, weights = NA, dim = 3)
```

## Arguments

- g:

  An `igraph` or a `tbl_graph` object

- pivots:

  Number of pivots to use for distance calculations

- weights:

  Optional edge weights to use for distance calculations

- dim:

  Desired number of dimensions for the layout

## Value

A matrix of coordinates for the pMDS layout
