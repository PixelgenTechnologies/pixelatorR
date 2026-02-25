# PolarizationScores

Get and set polarization scores for a
[`CellGraphAssay`](CellGraphAssay-class.md),
[`CellGraphAssay5`](CellGraphAssay5-class.md) or a `Seurat` object

## Usage

``` r
PolarizationScores(object, ...)

PolarizationScores(object, ...) <- value

# S3 method for class 'MPXAssay'
PolarizationScores(object, add_marker_counts = FALSE, ...)

# S3 method for class 'MPXAssay'
PolarizationScores(object, ...) <- value

# S3 method for class 'Seurat'
PolarizationScores(
  object,
  assay = NULL,
  meta_data_columns = NULL,
  add_marker_counts = FALSE,
  ...
)

# S3 method for class 'Seurat'
PolarizationScores(object, assay = NULL, ...) <- value
```

## Arguments

- object:

  An object with polarization scores

- ...:

  Not implemented

- value:

  A `tbl_df` with polarization scores

- add_marker_counts:

  A logical value indicating whether to add marker counts to the
  polarization score table.

- assay:

  Name of a `CellGraphAssay`

- meta_data_columns:

  A character vector with meta.data column names. This option can be
  useful to join meta.data columns with the polarization score table.

## Value

`PolarizationScores`: Polarization scores

`PolarizationScores<-`: An object with polarization scores updated

## See also

Other spatial metrics:
[`ColocalizationScores()`](ColocalizationScores.md),
[`Edgelists()`](Edgelists.md), [`ProximityScores()`](ProximityScores.md)
