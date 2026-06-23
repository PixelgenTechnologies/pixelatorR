# Compute approximate node saturation

Computes an approximate node saturation for each component in the
edgelist using SQL queries on a DuckDB connection. Node saturation
measures how many unique nodes (UMIs) have been observed relative to the
total number of reads.

## Usage

``` r
approximate_node_saturation(
  db,
  components = NULL,
  table_name = NULL,
  node_reads_multiplier = 2
)
```

## Arguments

- db:

  A `PixelDB` object with an active DuckDB connection.

- components:

  Optional character vector of component names to filter the edgelist.
  If `NULL` (default), all components are included.

- table_name:

  Optional name for the computed remote table. If provided, the result
  is materialized as a table in the DuckDB connection.

- node_reads_multiplier:

  Numeric; multiplier for the read count when calculating node
  saturation. Default is `2`, reflecting that each edge contributes to
  two nodes.

## Value

A lazy `tbl` (DuckDB query) with columns:

- component:

  The component (cell) identifier.

- nodes:

  The number of unique nodes.

- read_count:

  The total read count for the component.

- node_saturation:

  The node saturation value (0 to 1).

## Details

Node saturation is calculated as: \$\$S\_{node} = 1 - \frac{N}{R}\$\$
where \\N\\ is the number of unique nodes and \\R\\ is the total read
count (counted twice per edge, once for each endpoint).

This function creates a temporary view called `nodes` in the DuckDB
connection that aggregates node statistics (read count and degree) from
the edgelist.

Higher saturation values indicate that additional sequencing reads are
unlikely to reveal new nodes, suggesting the sample has been sequenced
to near-saturation.

## See also

[`approximate_edge_saturation`](approximate_edge_saturation.md) for
edge-based saturation.
[`approximate_saturation_curve`](approximate_saturation_curve.md) for
computing saturation at multiple downsampling fractions.
