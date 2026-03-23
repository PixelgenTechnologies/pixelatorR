# Calculate dispersion of components

Calculate the dispersion of numeric data in a Seurat object or matrix.
The dispersion is calculated as a summary statistic of the distribution
of values (such as counts) across components. The dispersion can be
calculated using the Gini coefficient or the Tau statistic. The Gini and
Tau coefficients are measures of inequality in a distribution, ranging
from 0 to 1, where 0 signifies a perfectly equal distribution and 1
signifies perfectly unequal distribution. In terms of marker counts
across different markers for a component, a low dispersion indicates
that the marker counts are at similar levels across all markers, while a
high dispersion indicates that the marker counts are distributed mostly
to one or a few markers. The dispersion can for example be used as a
quantitative metric to assess whether a component's counts are dispersed
between markers at an expected level. Both a high and low dispersion can
indicate that the component is some kind of an artifact, rather than a
real cell.

## Usage

``` r
CalculateDispersion(object, ...)

# S3 method for class 'matrix'
CalculateDispersion(object, method = c("gini", "tau"), margin = 2, ...)

# S3 method for class 'data.frame'
CalculateDispersion(object, method = c("gini", "tau"), margin = 2, ...)

# S3 method for class 'Matrix'
CalculateDispersion(object, method = c("gini", "tau"), margin = 2, ...)

# S3 method for class 'Seurat'
CalculateDispersion(
  object,
  method = c("gini", "tau"),
  margin = 2,
  assay = NULL,
  layer = NULL,
  metadata_name = NULL,
  ...
)
```

## Arguments

- object:

  A `Seurat` object or a count matrix.

- ...:

  Additional arguments. Currently not used.

- method:

  The method to use for calculating the dispersion. Can be "gini" or
  "tau". Default is "gini".

- margin:

  The margin of the object to apply the dispersion calculation on. 1
  indicates rows (typically markers), 2 indicates columns (typically
  cells). Default is 2. If `object` is a `Seurat` object and the margin
  is 2, a `Seurat` object is returned with the dispersion added as a
  metadata column.

- assay:

  A character with the name of the assay to use.

- layer:

  A character with the name of the layer to use.

- metadata_name:

  A character with the name of the metadata column to add. Default is
  `"dispersion_[method]"`.

## Value

A `Seurat` object with the dispersion added as a metadata column if
`margin = 2`, otherwise a numeric vector with the dispersion values for
each row or column of the input object.
