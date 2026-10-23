# Compute marginal discovery rate and saturation curve

Computes the instantaneous slope of the node and edge discovery curves
(\\dN / dR\\ and \\dE / dR\\) and the corresponding marginal
(Good-Turing) saturation across a vector of downsampling fractions
`fracs`.

## Usage

``` r
marginal_saturation_curve(
  db,
  fracs = seq(0.04, 1, by = 0.04),
  components = NULL,
  min_degree = 1L,
  verbose = TRUE
)
```

## Arguments

- db:

  A `PixelDB` object with an active DuckDB connection.

- fracs:

  A numeric vector of downsampling fractions in \\(0, 1\]\\. Default is
  `seq(0.04, 1.0, by = 0.04)`.

- components:

  Optional character vector of component names to include. If `NULL`
  (default), all components are included.

- min_degree:

  Integer (`1L` or `2L`); minimum dynamic node degree at fraction `p`
  for a node to be counted. Default is `1L` (all nodes). Set to `2L` to
  track dynamic degree \>= 2 core nodes.

- verbose:

  Logical; if `TRUE` (default), print progress messages.

## Value

A tibble with columns:

- component:

  The component (cell) identifier.

- read_count:

  Estimated number of reads at fraction `p`.

- p:

  The downsampling fraction.

- nodes:

  Expected number of unique nodes (with dynamic degree \>= `min_degree`)
  at fraction `p`.

- edges:

  Expected number of unique edges at fraction `p`.

- degree:

  Estimated average node degree (\\2 \times \text{edges} /
  \text{nodes}\\).

- node_discovery_rate:

  Instantaneous slope \\dN / dR\\ at fraction `p` (new nodes per
  additional read).

- edge_discovery_rate:

  Instantaneous slope \\dE / dR\\ at fraction `p` (new edges per
  additional read).

- degree_rate:

  Instantaneous slope \\d\bar{d} / dR\\ at fraction `p` (change in
  average node degree per additional read).

- marginal_node_saturation:

  Marginal node saturation at fraction `p` (\\1 - (dN / dR) / 2\\).

- marginal_edge_saturation:

  Marginal (Good-Turing) edge saturation at fraction `p`.

## Details

### Marginal discovery rate across downsampling fractions

At any downsampling fraction \\p \in (0, 1\]\\, the expected number of
retained edges is \\\mathbb{E}\[E(p)\] = \sum\_{e=1}^{E} (1 - (1 -
p)^{r_e})\\. Differentiating with respect to the downsampled read count
\\R(p) = p R\\ gives the exact marginal edge discovery rate (new edges
per additional read): \$\$\text{edge\\discovery\\rate}(p) =
\frac{d\\\mathbb{E}\[E(p)\]}{dR(p)} = \frac{1}{R} \sum\_{e=1}^{E} r_e
(1 - p)^{r_e - 1}\$\$

For nodes:

- When `min_degree = 1L`, \\\mathbb{E}\[N(p)\] = \sum\_{v=1}^{N} (1 -
  (1 - p)^{r_v})\\, so the exact marginal node discovery rate is:
  \$\$\text{node\\discovery\\rate}(p) =
  \frac{d\\\mathbb{E}\[N(p)\]}{dR(p)} = \frac{1}{R} \sum\_{v=1}^{N} r_v
  (1 - p)^{r_v - 1}\$\$

- When `min_degree = 2L`, a node has dynamic degree \\\ge 2\\ at
  fraction \\p\\ with exact probability: \$\$P(\text{deg}\_v(p) \ge 2) =
  1 + (d_v - 1)(1 - p)^{r_v} - \sum\_{e \ni v} (1 - p)^{r_v - r_e}\$\$
  Differentiating with respect to \\R(p) = p R\\ yields the exact
  marginal discovery rate of dynamic degree \\\ge 2\\ core nodes:
  \$\$\text{node\\discovery\\rate}(p) = \frac{1}{R} \sum\_{v:\\ d_v \ge
  2} \left\[ \sum\_{e \ni v} (r_v - r_e)(1 - p)^{r_v - r_e - 1} -
  (d_v - 1) r_v (1 - p)^{r_v - 1} \right\]\$\$

## See also

[`marginal_saturation`](marginal_saturation.md) for computing full-depth
marginal discovery rates and saturation as a lazy DuckDB table.
[`approximate_saturation_curve`](approximate_saturation_curve.md) for
cumulative saturation curves.
