# Perform spatial smoothing of a count matrix using the adjacency matrix of the graph

Perform spatial smoothing of a count matrix using the adjacency matrix
of the graph

## Usage

``` r
spatial_smoothing(A, x, iter = 5L)
```

## Arguments

- A:

  The adjacency matrix of the graph (sparse matrix of class
  `dgCMatrix`).

- x:

  The count matrix to be smoothed (sparse matrix of class `dgCMatrix`).

- iter:

  The number of iterations to perform the smoothing. Default is 5.

## Value

A smoothed count matrix of the same dimensions as `x`.
