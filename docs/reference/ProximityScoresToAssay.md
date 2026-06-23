# Convert proximity score table to an Assay or Assay5

Convert proximity score table to an Assay or Assay5

## Usage

``` r
ProximityScoresToAssay(object, ...)

# S3 method for class 'tbl_lazy'
ProximityScoresToAssay(
  object,
  values_from = "log2_ratio",
  separator = ":",
  ...
)

# S3 method for class 'data.frame'
ProximityScoresToAssay(
  object,
  values_from = "log2_ratio",
  separator = ":",
  ...
)

# S3 method for class 'PNAAssay'
ProximityScoresToAssay(
  object,
  values_from = "log2_ratio",
  separator = ":",
  lazy = FALSE,
  ...
)

# S3 method for class 'PNAAssay5'
ProximityScoresToAssay(
  object,
  values_from = "log2_ratio",
  separator = ":",
  lazy = FALSE,
  ...
)

# S3 method for class 'Seurat'
ProximityScoresToAssay(
  object,
  assay = NULL,
  new_assay = NULL,
  values_from = "log2_ratio",
  separator = ":",
  lazy = FALSE,
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
  pick values from. Default is "log2_ratio".

- separator:

  A character to separate marker names in the row names of the output.
  Default is ":". Must be a single character and must not appear in any
  marker name.

- lazy:

  Whether to look for proximity scores in the PXL file instead of the
  [`Assay`](https://satijalab.github.io/seurat-object/reference/Assay-class.html)
  /
  [`Assay5`](https://satijalab.github.io/seurat-object/reference/Assay5-class.html)
  object.

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
must contain the following columns:

- "marker_1": the name of the first marker in the pair

- "marker_2": the name of the second marker in the pair

- "component": the name of the component

- a column with the proximity scores to be used as values in the wide
  format (defined by the `values_from` parameter)

The wide format is an array-like object with dimensions pair x
components, where pair is defined as the combination of "marker_1" and
"marker_2" separated by the `separator`, and where each element in the
array is filled with a value for a selected spatial metric.

Note that observations that are missing from the proximity score table
are replaced with 0's. Proximity scores can also be 0 (no deviation from
random expectations), and it will not be possible to distinguish between
these two cases in the output.

Different outputs are returned depending on the input object type:

- `tibble/data.frame`: returns a `dgCMatrix` with marker pairs in rows
  and components in columns. The components are not ordered in any
  particular way and should therefore be ordered before placing them in
  e.g. a Seurat object.

- `PNAAssay/PNAAssay5`: returns an `Assay` or `Assay5` with marker pairs
  in rows and components in columns. Columns are ordered according to
  the order of the components in the `PNAAssay/PNAAssay5`.

- `Seurat` object: returns the `Seurat` object with a new `Assay` or
  `Assay5` with marker pairs in rows and components in columns. Columns
  are ordered according to the order of the components in the `Seurat`
  object.

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
