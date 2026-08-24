# Approximate nearest neighbors using Annoy

This function takes a matrix, and returns approximate Euclidean nearest
neighbors and distances of row items given the number of trees (n_trees)
and number of nearest neighbors (n_nn).

## Usage

``` r
FindAnnoyNeighbors(
  x,
  cells = NULL,
  n_trees = 50L,
  n_nn = 10L,
  search_k = NULL,
  annoy_alg = c("euclidean", "angular", "manhattan", "hamming")
)
```

## Arguments

- x:

  A numeric matrix with data to find nearest neighbors. Rows are cells,
  and columns are features.

- cells:

  A character vector with cell names to find nearest neighbors for. If
  NULL, all cells are used.

- n_trees:

  An integer with the number of trees to build in the Annoy index.

- n_nn:

  An integer with the number of nearest neighbors to find.

- search_k:

  An integer with the number of nodes to search in the Annoy index.
  Default is `n_trees * n_nn`.

- annoy_alg:

  An character specifying which distance algorithm to use. Default is
  [`AnnoyEuclidean`](https://rdrr.io/pkg/RcppAnnoy/man/AnnoyIndex.html)
  (`"euclidean"`). Available options are `"euclidean"`, `"angular"`,
  `"manhattan"`, and `"hamming"`.

## Value

A tibble with the following columns:

- id: The cell name

- index: The index of the cell in the Annoy index

- item: The index of the nearest neighbor

- distance: The distance (Euclidean by default) to the nearest neighbor

- nn: The rank of the nearest neighbor

- neighbor: The name of the nearest neighbor

## Examples

``` r
x <- matrix(rnorm(1000), ncol = 10)
FindAnnoyNeighbors(x, n_trees = 50, n_nn = 10)
#> # A tibble: 1,000 × 6
#>    id    index  item distance    nn neighbor
#>    <chr> <dbl> <int>    <dbl> <int> <chr>   
#>  1 1         0     0     0        1 1       
#>  2 1         0    12     2.51     2 13      
#>  3 1         0    44     2.77     3 45      
#>  4 1         0    29     2.77     4 30      
#>  5 1         0     8     2.82     5 9       
#>  6 1         0    47     2.99     6 48      
#>  7 1         0    22     3.11     7 23      
#>  8 1         0    24     3.16     8 25      
#>  9 1         0    40     3.26     9 41      
#> 10 1         0    31     3.35    10 32      
#> # ℹ 990 more rows
```
