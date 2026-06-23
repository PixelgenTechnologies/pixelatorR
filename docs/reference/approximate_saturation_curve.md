# Compute approximate saturation curve

Computes an approximate saturation curve for nodes and edges in an
edgelist by estimating the expected number of retained elements at
various downsampling fractions. This provides insight into how
sequencing depth affects saturation without requiring actual
downsampling.

## Usage

``` r
approximate_saturation_curve(
  db,
  fracs = seq(0.04, 0.96, by = 0.04),
  components = NULL,
  node_reads_multiplier = 2,
  verbose = TRUE
)
```

## Arguments

- db:

  A `PixelDB` object with an active DuckDB connection.

- fracs:

  A numeric vector of fractions to downsample the edgelist. Values must
  be strictly between 0 and 1. Default is `seq(0.04, 0.96, by = 0.04)`.

- components:

  Optional character vector of component names to compute saturation
  for. If `NULL` (default), all components are included.

- node_reads_multiplier:

  Numeric; multiplier for the read count when calculating node
  saturation. Default is `2`, reflecting that each edge contributes to
  two nodes. Optionally, set to `1` to use the same read count for nodes
  and edges, which may be more appropriate in certain contexts.

- verbose:

  Logical; if `TRUE` (default), print progress messages.

## Value

A tibble with the following columns:

- component:

  The component (cell) identifier.

- read_count:

  The estimated number of reads after downsampling.

- p:

  The downsampling fraction.

- nodes:

  The estimated number of nodes retained.

- edges:

  The estimated number of edges retained.

- degree:

  The estimated average node degree.

- node_saturation:

  The node saturation at the given fraction.

- edge_saturation:

  The edge saturation at the given fraction.

## Details

### Downsampling model

The downsampling is performed probabilistically using expected values
rather than actual random sampling. For each edge with a given
`read_count`, the probability of the edge being retained at fraction `p`
is: \$\$P(\text{edge retained}) = 1 - (1 - p)^{\text{read\\count}}\$\$

For nodes, the probability of retention depends on the node's degree. A
node is retained if at least one of its incident edges is retained:
\$\$P(\text{node retained}) = 1 - (1 - p\_{edges})^{\text{degree}}\$\$
where \\p\_{edges}\\ is the fraction of edges retained.

### Saturation calculation

The saturation values are calculated using: \$\$S = 1 -
\frac{N\_{downsampled}}{R\_{downsampled}}\$\$

where \\N\_{downsampled}\\ is the number of nodes or edges in the
downsampled graph, and \\R\_{downsampled}\\ is the number of reads
supporting those elements. Since each read supports two nodes (one at
each endpoint), \\R\_{downsampled, nodes} = 2 \times R\_{downsampled,
edges}\\.

### Implementation notes

This function requires that
[`approximate_node_saturation`](approximate_node_saturation.md) has been
called first (either directly or through this function) to create the
`nodes` temporary view in the DuckDB connection.

## See also

[`approximate_node_saturation`](approximate_node_saturation.md) and
[`approximate_edge_saturation`](approximate_edge_saturation.md) for
computing saturation at full depth. [`lcc_curve`](lcc_curve.md) for
computing largest connected component sizes at different downsampling
fractions.

## Examples

``` r
if (FALSE) { # \dontrun{
library(ggplot2)
library(dplyr)

pxl_file <- minimal_pna_pxl_file()
db <- PixelDB$new(pxl_file)

sat_curve <- approximate_saturation_curve(db, fracs = seq(0.1, 0.9, by = 0.1))

# Plot node saturation curve
sat_curve %>%
  tidyr::pivot_longer(
    cols = c(node_saturation, edge_saturation),
    names_to = "type",
    values_to = "saturation"
  ) %>%
  ggplot(aes(read_count, saturation, color = component)) +
  geom_line() +
  geom_point() +
  labs(x = "Reads", y = "Saturation") +
  theme_minimal() +
  facet_grid(~type) +
  scale_y_continuous(limits = c(0, 1))

db$close()
} # }
```
