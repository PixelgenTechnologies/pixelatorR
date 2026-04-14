# Compute transition probabilities

Create a transition probability matrix from an adjacency matrix.

## Usage

``` r
compute_transition_probabilities(A, k = 1, remove_self_loops = FALSE)
```

## Arguments

- A:

  An adjacency matrix

- k:

  **\[experimental\]** An integer value specifying the maximum steps
  from each node to expand the neighborhood. E.g. with `k = 2`, each
  node neighborhood will include the direct neighbors and nodes two
  steps away. With `k = 1` (default), only the direct neighbors will be
  used.

- remove_self_loops:

  Whether to remove self-loops from the transition probability matrix.

## Value

A matrix of transition probabilities. The transition probability \\P^k(u
\rightarrow v)\\ for a k-step walk is found in row \\u\\ and column
\\v\\ of the transition probability matrix. Row \\v\\ and column \\u\\
gives the reversed transition probability \\P^k(v \rightarrow u)\\ for
the k-step walk.

## Examples

``` r
library(igraph)
#> Warning: package 'igraph' was built under R version 4.5.3
#> 
#> Attaching package: 'igraph'
#> The following object is masked from 'package:tidygraph':
#> 
#>     groups
#> The following objects are masked from 'package:dplyr':
#> 
#>     as_data_frame, groups, union
#> The following objects are masked from 'package:stats':
#> 
#>     decompose, spectrum
#> The following object is masked from 'package:base':
#> 
#>     union
g <- make_lattice(c(2, 3))
A <- as_adjacency_matrix(g)

# Transition probabilities
P_out <- compute_transition_probabilities(A)
P_out
#> 6 x 6 sparse Matrix of class "dgCMatrix"
#>                                                                 
#> [1,] .         0.5000000 0.5000000 .         .         .        
#> [2,] 0.5000000 .         .         0.5000000 .         .        
#> [3,] 0.3333333 .         .         0.3333333 0.3333333 .        
#> [4,] .         0.3333333 0.3333333 .         .         0.3333333
#> [5,] .         .         0.5000000 .         .         0.5000000
#> [6,] .         .         .         0.5000000 0.5000000 .        
```
