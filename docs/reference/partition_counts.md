# Compute partition counts

Compute partition counts

## Usage

``` r
partition_counts(cg, partition = NULL, partition_column = NULL)
```

## Arguments

- cg:

  A CellGraph object containing the cell graph and count data.

- partition:

  A character or factor vector indicating the partition of nodes into
  groups. Either `partition` or `partition_column` must be provided.

- partition_column:

  A string indicating the name of the vertex attribute in the cell graph
  that contains the partition information.

## Value

A matrix of counts for each partition group, where rows correspond to
partition groups and columns correspond to features (e.g., proteins).
