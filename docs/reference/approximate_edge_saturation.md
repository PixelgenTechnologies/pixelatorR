# Compute approximate edge saturation

Computes an approximate edge saturation for each component in the
edgelist using SQL queries on a DuckDB connection. Edge saturation
measures how many unique edges have been observed relative to the total
number of reads.

## Usage

``` r
approximate_edge_saturation(db, components = NULL, table_name = NULL)
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

## Value

A lazy `tbl` (DuckDB query) with columns:

- component:

  The component (cell) identifier.

- edges:

  The number of unique edges.

- read_count:

  The total read count for the component.

- edge_saturation:

  The edge saturation value (0 to 1).

## Details

Edge saturation is calculated as: \$\$S\_{edge} = 1 - \frac{E}{R}\$\$
where \\E\\ is the number of unique edges and \\R\\ is the total read
count.

Higher saturation values indicate that additional sequencing reads are
unlikely to reveal new edges, suggesting the sample has been sequenced
to near-saturation.

## See also

[`approximate_node_saturation`](approximate_node_saturation.md) for
node-based saturation.
[`approximate_saturation_curve`](approximate_saturation_curve.md) for
computing saturation at multiple downsampling fractions.
