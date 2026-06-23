# Calculate antibody counts per A-node

**\[deprecated\]**

Computes and returns a data frame of antibody counts per node (vertex)
of the A node graph given a component edge list as input. The parameter
`k` allows to include neighbors (of each node) when computing the
counts. `k` defines the number of levels when searching neighbors.

If `k > 0`, [`edgelist_to_simple_Anode_graph`](graph-conversion.md) is
used to calculate an A-node projected graph.

## Usage

``` r
node_markers_counts(component_edge_list, k = 0)
```

## Arguments

- component_edge_list:

  An object of class `tbl_df`

- k:

  Number of neighbors to include

## Value

A matrix with node marker counts
