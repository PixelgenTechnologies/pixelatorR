# Scale layout coordinates

This function scales the layout coordinates to a given percentile of the
maximum radius.

## Usage

``` r
scale_layout(layout, percentile = 0.95)
```

## Arguments

- layout:

  A tibble with x, y and optionally z coordinates.

- percentile:

  A numeric value between 0 and 1 indicating the percentile of the
  maximum radius to scale the layout to.

## Value

A tibble with the scaled coordinates.
