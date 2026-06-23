# Add node colors to a `CellGraph`

**\[deprecated\]**

Adds node colors based on the expression of a single or multiple
markers. If multiple markers are provided, their values are either
multiplied or summed depending on the `mode`.

## Usage

``` r
color_by_marker(
  cg,
  markers,
  smooth_counts = TRUE,
  palette = "viridis",
  rev_pal = FALSE,
  mode = c("product", "sum"),
  normalize = TRUE,
  trim_quantiles = c(0, 1),
  nNodes = NULL
)
```

## Arguments

- cg:

  A [`CellGraph`](CellGraph-class.md) object

- markers:

  A character vector with valid marker IDs

- smooth_counts:

  Applies neighborhood smoothing of marker counts. Each node marker
  counts is combined with the sum of marker counts in a neighborhood of
  length 1.

- palette:

  A valid color palette from `RColorBrewer` or `viridis` or a vector of
  colors.

- rev_pal:

  Should the palette be reversed?

- mode:

  Only used if more than 1 marker is provided. With `mode = "product"`,
  the marker counts are multiplied and with `mode = "sum"`, the marker
  counts are summed. This step is executed after smoothing.

- normalize:

  Apply normalization to the marker counts after smoothing and the
  multiplication/summation step. The values are first divided by the
  node degree followed by log-transformation with `log1p`.

- trim_quantiles:

  A numeric vector of length 2 specifying quantiles to trim. This can be
  useful to reduce the influence of outliers. The specified quantiles
  are calculated from the marker counts values prior to normalization.
  (see section below for details)

- nNodes:

  Number of nodes to keep in the graph

## Trim quantile

Sometimes, a few outliers can dominate which makes the colors difficult
to interpret. With `trim_quantiles`, we can remove such outliers. By
setting the upper bound to 0.99, we are calculating the 99th quantile of
the marker count values and any values above this threshold are replaced
with the value for the 99th quantile. The same thing applies to the
lower bound but the other way around.
