# Expand an adjacency matrix to include higher-order neighborhoods

The node neighborhoods are expanded by multiplying the adjacency matrix
`k` times. To make sure that all neighbors up to `k` steps away are
included, self-loops are added by setting the diagonal of `A` to 1.
These self-loops are removed from the resulting matrix after expansion.
Note that the self-loops affect how the transition probabilities are
calculated. When `use_weights = FALSE`, the edge values of the result
from the matrix multiplication corresponds to the number of k-step
random walks between any pair of nodes. For the expanded adjacency
matrix, we replace these edge values with 1 so that the adjacency matrix
corresponds to an unweighted graph with expanded neighborhoods.

## Usage

``` r
expand_adjacency_matrix(A, k = 1L, use_weights = FALSE, min_weight = 0)
```

## Arguments

- A:

  An adjacency matrix. Can be a sparse matrix (`dgCMatrix`)

- k:

  An integer specifying the neighborhood size. Each node will be
  connected to nodes up to `k` steps away.

- use_weights:

  If `TRUE`, the transpose of the stochastic matrix will be expanded.
  The resulting edge weights will correspond to the incoming transition
  probabilities for a `k`-step random walk.

- min_weight:

  A numeric value between 0 and 1 specifying the minimum weight
  threshold for edges in the expanded adjacency matrix. Only applicable
  when `use_weights = TRUE`. Edges with weights below this threshold
  will be set to zero, effectively removing them from the graph. This
  can help to reduce noise and focus on stronger connections in the
  expanded neighborhood.

## Value

An expanded adjacency matrix
