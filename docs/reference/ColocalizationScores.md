# ColocalizationScores

Get and set colocalization scores for a
[`CellGraphAssay`](CellGraphAssay-class.md),
[`CellGraphAssay5`](CellGraphAssay5-class.md) or a `Seurat` object

## Usage

``` r
ColocalizationScores(object, ...)

ColocalizationScores(object, ...) <- value

# S3 method for class 'MPXAssay'
ColocalizationScores(object, add_marker_counts = FALSE, ...)

# S3 method for class 'MPXAssay'
ColocalizationScores(object, ...) <- value

# S3 method for class 'Seurat'
ColocalizationScores(
  object,
  assay = NULL,
  meta_data_columns = NULL,
  add_marker_counts = FALSE,
  ...
)

# S3 method for class 'Seurat'
ColocalizationScores(object, assay = NULL, ...) <- value
```

## Arguments

- object:

  An object with polarization scores

- ...:

  Not implemented

- value:

  A `tbl_df` with colocalization scores

- add_marker_counts:

  A logical value indicating whether to add marker counts to the
  colocalization score table.

- assay:

  Name of a `CellGraphAssay`

- meta_data_columns:

  A character vector with meta.data column names. This option can be
  useful to join meta.data columns with the colocalization score table.

## Value

`ColocalizationScores`: Colocalization scores

`ColocalizationScores<-`: An object with colocalization scores updated

## See also

Other spatial metrics: [`Edgelists()`](Edgelists.md),
[`PolarizationScores()`](PolarizationScores.md),
[`ProximityScores()`](ProximityScores.md)
