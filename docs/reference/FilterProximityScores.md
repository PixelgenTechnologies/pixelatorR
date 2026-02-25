# Filter proximity scores

At least one of the thresholds must be set. The function will filter out
protein pairs that do not meet the specified thresholds.

All filters can be used together or separately, but it's recommended to
only use the `background_threshold_pct` filter. This filter is applied
to the fraction of UMI counts which does not depend on the cell size
(total number of molecules).

The filters are applied in the following order:

1.  Remove protein pairs with a minimum UMI count fraction below
    `background_threshold_pct`

2.  Remove protein pairs with a minimum UMI count below
    `background_threshold_count`

3.  Remove protein pairs detected in fewer than `min_cells_count` cells

## Usage

``` r
FilterProximityScores(object, ...)

# S3 method for class 'tbl_lazy'
FilterProximityScores(
  object,
  background_threshold_pct = NULL,
  background_threshold_count = NULL,
  min_cells_count = NULL,
  name = "filtered_proximity",
  ...
)

# S3 method for class 'data.frame'
FilterProximityScores(
  object,
  background_threshold_pct = NULL,
  background_threshold_count = NULL,
  min_cells_count = NULL,
  ...
)
```

## Arguments

- object:

  An `tbl_df` or `tbl_lazy` object with proximity scores.

- ...:

  Additional arguments. Currently not used.

- background_threshold_pct:

  The background abundance level given as a fraction total UMI counts.
  For example, 0.05 means that the background is set at 5% of the total
  UMI counts for each cell. A protein pair must have at least this
  fraction for both proteins to be kept in the table.

- background_threshold_count:

  The background abundance level given as a UMI count. For example, 30
  means that the background is set at 30 UMI counts. A protein pair must
  have at least this UMI count for both proteins to be kept in the
  table. Note that this puts a hard cutoff on the UMI counts, which may
  not be desired if the cells have very different total UMI counts.

- min_cells_count:

  The minimum number of cells in which a protein pair must be detected
  to be kept in the table. This is useful to remove protein pairs that
  are detected in very few cells, which may be due to noise or low
  expression levels.

- name:

  Name of the temporary table in the database containing the filtered
  proximity scores.

## Value

A `tbl_df` or `tbl_lazy` with filtered proximity scores.

## Examples

``` r
library(pixelatorR)

pxl_file <- minimal_pna_pxl_file()
se <- ReadPNA_Seurat(pxl_file)
#> ✔ Created a <Seurat> object with 5 cells and 158 targeted surface proteins
proximity_table <- ProximityScores(se, add_marker_proportions = TRUE)
#> ! Setting `add_marker_counts = TRUE` which is required when `add_marker_proportions = TRUE`.

# Filter scores
proximity_table_filtered <- FilterProximityScores(
  proximity_table,
  background_threshold_pct = 0.001
)
#> Warning: Missing values are always removed in SQL aggregation functions.
#> Use `na.rm = TRUE` to silence this warning
#> This warning is displayed once every 8 hours.

# Rows kept
pct_rows_kept <- round(nrow(proximity_table_filtered) / nrow(proximity_table) * 100, digits = 2)
glue::glue("Fraction of rows kept: {pct_rows_kept}%")
#> Fraction of rows kept: 11.13%
```
