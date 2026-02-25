# Compute approximate saturation curve

This function computes an approximate saturation curve for nodes and
edges in an edgelist. The edgelist is downsampled to various fractions
of the total number of reads, and the number of remaining nodes and
edges are calculated for each fraction. The saturation values are then
calculated relative to the theoretical maximum number of nodes and edges
using the Chao1 estimator. The theoretical maximum is computed for each
component in the full edgelist. Note that Chao1 will give a lower bound
for the theoretical maximum, hence the saturation values are likely
overestimated. The estimate will be more robust if the sample was
sequenced at a high depth.

## Usage

``` r
approximate_saturation_curve(
  db,
  fracs = seq(0.04, 0.96, by = 0.04),
  detailed = FALSE,
  verbose = TRUE
)
```

## Arguments

- db:

  A `PixelDB` object.

- fracs:

  A numeric vector of fractions to downsample the edgelist.

- detailed:

  If `TRUE`, the function will return ...

- verbose:

  Print messages

## Value

A tibble with the saturation values
