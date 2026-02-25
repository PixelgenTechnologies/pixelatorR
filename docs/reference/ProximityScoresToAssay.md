# Convert proximity score table to an Assay or Assay5

Convert proximity score table to an Assay or Assay5

## Usage

``` r
ProximityScoresToAssay(object, ...)

# S3 method for class 'tbl_lazy'
ProximityScoresToAssay(object, values_from = "join_count_z", ...)

# S3 method for class 'data.frame'
ProximityScoresToAssay(
  object,
  values_from = "join_count_z",
  missing_obs = NA_real_,
  return_sparse = TRUE,
  ...
)

# S3 method for class 'PNAAssay'
ProximityScoresToAssay(
  object,
  values_from = "join_count_z",
  missing_obs = NA_real_,
  ...
)

# S3 method for class 'PNAAssay5'
ProximityScoresToAssay(
  object,
  values_from = "join_count_z",
  missing_obs = NA_real_,
  ...
)

# S3 method for class 'Seurat'
ProximityScoresToAssay(
  object,
  assay = NULL,
  new_assay = NULL,
  values_from = "join_count_z",
  missing_obs = NA_real_,
  ...
)
```

## Arguments

- object:

  An object with proximity scores

- ...:

  Not yet implemented

- values_from:

  A single string defining what column in the proximity score table to
  pick values from.

- missing_obs:

  A numeric value or NA to replace missing observations with. Default is
  `NA_real_`.

- return_sparse:

  A logical specifying whether to return a sparse matrix (`dgCMatrix`).

- assay:

  Name of the [`PNAAssay`](PNAAssay-class.md) or
  [`PNAAssay5`](PNAAssay5-class.md) to pull proximity scores from

- new_assay:

  Name of the
  [`Assay`](https://satijalab.github.io/seurat-object/reference/Assay-class.html)
  or
  [`Assay5`](https://satijalab.github.io/seurat-object/reference/Assay5-class.html)
  to store the wide formatted spatial metric in

## Behavior

Takes an object with PNA proximity scores in long format and returns an
object with proximity scores in a wide format. The proximity score table
contains various spatial metrics along with p-values for each protein
pair and component.

The wide format is an array-like object with dimensions (markers_1 \*
marker_2) x components, where each cell is filled with a value for a
selected spatial metric.

Note that that observations that are missing from the proximity score
table can be replaced with 0's in the wide array by setting
`missing_obs_val = 0`, which might be required by various functions in
`Seurat`. However, this may not be the desired behavior as a value of 0
usually doesn't mean that the observation is missing.

Different outputs are returned depending on the input object type:

- `tibble/data.frame`: returns a `matrix` with marker pairs in rows and
  components in columns

- `PNAAssay/PNAAssay5`: returns an `Assay` or `Assay5` with marker pairs
  in rows and components in columns

- `Seurat` object: returns the `Seurat` object with a new `Assay` or
  `Assay5` with marker pairs in rows and components in columns

As many methods provided in Seurat operates on `Assay`/`Assay5` objects,
it can sometimes be convenient to make this conversion if you wish to
use these methods on spatial metrics in the proximity score table. For
instance, if we want to compute a UMAP on the proximity scores with
`RunUMAP`, we need the values to be formatted in an `Assay`/`Assay5`.
This also makes it possible to use various visualization functions such
as `VlnPlot` or `FeaturePlot` to show the distribution of proximity
scores.

## See also

Other Spatial metrics conversion methods:
[`ColocalizationScoresToAssay()`](ColocalizationScoresToAssay.md),
[`PolarizationScoresToAssay()`](PolarizationScoresToAssay.md)
