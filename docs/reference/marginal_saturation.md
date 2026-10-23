# Compute marginal discovery rate and saturation

Computes the marginal discovery rate (the instantaneous slope of the
discovery curve, i.e., expected new nodes and edges discovered per
additional sequencing read) and the corresponding marginal (Good-Turing)
saturation at full sequencing depth for each component in a `PixelDB`
object.

## Usage

``` r
marginal_saturation(db, components = NULL, table_name = NULL, min_degree = 1L)
```

## Arguments

- db:

  A `PixelDB` object with an active DuckDB connection.

- components:

  Optional character vector of component names to include. If `NULL`
  (default), all components are included.

- table_name:

  Optional name for materializing the result as a table in the DuckDB
  connection.

- min_degree:

  Integer (`1L` or `2L`); minimum node degree at the current depth for a
  node to be counted. Default is `1L` (all nodes). Set to `2L` to
  measure the discovery rate of dynamic degree \>= 2 core nodes.

## Value

A lazy `tbl` (DuckDB query) with columns:

- component:

  The component (cell) identifier.

- read_count:

  Total sequencing reads in the component.

- nodes:

  Number of unique nodes (with degree \>= `min_degree`).

- edges:

  Number of unique edges in the component.

- degree:

  Average node degree (\\2 \times \text{edges} / \text{nodes}\\).

- node_discovery_rate:

  Instantaneous slope \\dN / dR\\ (new nodes per additional read).

- edge_discovery_rate:

  Instantaneous slope \\dE / dR\\ (new edges per additional read).

- degree_rate:

  Instantaneous slope \\d\bar{d} / dR\\ (change in average node degree
  per additional read).

- marginal_node_saturation:

  Marginal node saturation (\\1 - (dN / dR) / 2\\).

- marginal_edge_saturation:

  Marginal (Good-Turing) edge saturation (\\1 - dE / dR\\).

## Details

### Marginal discovery rate (slope)

Let \\\mathbb{E}\[U(p)\]\\ be the expected number of unique elements
(nodes or edges) observed at downsampling fraction \\p \in (0, 1\]\\,
and let \\R(p) = p R\\ be the downsampled read count. The instantaneous
slope of the discovery curve with respect to sequencing reads is:
\$\$\frac{d\\\mathbb{E}\[U(p)\]}{dR(p)} = \frac{1}{R}
\frac{d\\\mathbb{E}\[U(p)\]}{dp}\$\$

Evaluated at full sequencing depth (\\p = 1\\):

- **Edges:** An edge with \\r_e\\ reads has retention probability \\1 -
  (1 - p)^{r_e}\\, whose derivative at \\p = 1\\ is \\1\\ if \\r_e = 1\\
  and \\0\\ if \\r_e \ge 2\\. Thus, the edge discovery rate at full
  depth is the fraction of reads supporting singleton edges (\\E_1\\):
  \$\$\text{edge\\discovery\\rate} =
  \left.\frac{d\\\mathbb{E}\[E\]}{dR}\right\|\_{p=1} = \frac{E_1}{R}\$\$

- **Nodes (`min_degree = 1`):** A node with total incident reads \\r_v\\
  has retention probability \\1 - (1 - p)^{r_v}\\, so at \\p = 1\\:
  \$\$\text{node\\discovery\\rate} =
  \left.\frac{d\\\mathbb{E}\[N\]}{dR}\right\|\_{p=1} = \frac{N_1}{R}\$\$
  where \\N_1\\ is the number of singleton nodes (\\r_v = 1\\).

- **Core nodes (`min_degree = 2`):** A node is in the dynamic 2-core at
  fraction \\p\\ if at least two of its incident edges are retained:
  \$\$P(\text{deg}\_v(p) \ge 2) = 1 + (d_v - 1)(1 - p)^{r_v} - \sum\_{e
  \ni v} (1 - p)^{r_v - r_e}\$\$ Taking the derivative with respect to
  \\p\\ at \\p = 1\\, the only terms with \\r_v - r_e = 1\\ occur when a
  degree-2 node (\\d_v = 2\\) has an incident singleton edge (\\r_e =
  1\\), so: \$\$\text{node\\discovery\\rate} =
  \left.\frac{d\\\mathbb{E}\[N\_{\text{deg} \ge 2}\]}{dR}\right\|\_{p=1}
  = \frac{\sum\_{v:\\ d_v = 2} s_v}{R}\$\$ where \\s_v \in \\0, 1, 2\\\\
  is the number of incident singleton edges on degree-2 node \\v\\.

Multiplying `node_discovery_rate` or `edge_discovery_rate` by `10000`
gives the expected number of new unique nodes or edges gained per 10,000
additional sequencing reads.

### Marginal (Good-Turing) saturation

The marginal saturation is the probability that the next observation
lands on an already-observed element rather than discovering a new one
(since each sequencing read observes 2 node endpoints,
`node_discovery_rate` is divided by `2`): \$\$S\_{\text{marginal, edge}}
= 1 - \text{edge\\discovery\\rate}\$\$ \$\$S\_{\text{marginal, node}} =
1 - \frac{\text{node\\discovery\\rate}}{2}\$\$

The instantaneous rate of change of average node degree (\\\bar{d} = 2E
/ N\\) per additional read is given by the quotient rule:
\$\$\text{degree\\rate} = \frac{d\bar{d}}{dR} =
\frac{2\\\text{edge\\discovery\\rate} -
\bar{d}\\\text{node\\discovery\\rate}}{N}\$\$

## See also

[`marginal_saturation_curve`](marginal_saturation_curve.md) for
computing marginal discovery rates and saturation across downsampling
fractions.
[`approximate_node_saturation`](approximate_node_saturation.md) and
[`approximate_edge_saturation`](approximate_edge_saturation.md) for
cumulative saturation.
