# Fetch graph layout coordinates

Fetches a stored 3D layout from [`CellGraph`](CellGraph-class.md)
objects and optionally joins node-level variables retrieved with
[`FetchData`](https://satijalab.github.io/seurat-object/reference/FetchData.html).

## Usage

``` r
FetchLayoutData(object, ...)

# S3 method for class 'CellGraph'
FetchLayoutData(
  object,
  layout_method = "wpmds_3d",
  vars = NULL,
  add_marker = FALSE,
  layer = NULL,
  ...
)

# S3 method for class 'CellGraphList'
FetchLayoutData(
  object,
  layout_method = "wpmds_3d",
  vars = NULL,
  cells = NULL,
  add_marker = FALSE,
  layer = NULL,
  ...
)

# S3 method for class 'PNAAssay'
FetchLayoutData(
  object,
  layout_method = "wpmds_3d",
  vars = NULL,
  cells = NULL,
  add_marker = FALSE,
  layer = NULL,
  ...
)

# S3 method for class 'PNAAssay5'
FetchLayoutData(
  object,
  layout_method = "wpmds_3d",
  vars = NULL,
  cells = NULL,
  add_marker = FALSE,
  layer = NULL,
  ...
)

# S3 method for class 'Seurat'
FetchLayoutData(
  object,
  layout_method = "wpmds_3d",
  vars = NULL,
  cells = NULL,
  assay = NULL,
  add_marker = FALSE,
  layer = NULL,
  ...
)
```

## Arguments

- object:

  An object

- ...:

  Additional parameters passed to other methods

- layout_method:

  Name of a stored layout, typically one computed with
  [`ComputeLayout`](ComputeLayout.md) or loaded with
  [`LoadCellGraphs`](LoadCellGraphs.md). Default is `"wpmds_3d"`. The
  layout must contain `x`, `y`, and `z` coordinates.

- vars:

  Optional character vector of node-level variables to fetch with
  [`FetchData`](https://satijalab.github.io/seurat-object/reference/FetchData.html)
  (markers, metadata columns, graph vertex attributes, or reduction
  embeddings). A variable that is present on some cells and missing on
  others is filled with `NA` and a warning names those cells. A variable
  missing from every cell aborts.

- add_marker:

  If `TRUE`, add a `marker` column with the marker label of each node.
  Labels are read from the counts matrix. Nodes with no count are `NA`.

- layer:

  Name of a node matrix layer passed to
  [`FetchData`](https://satijalab.github.io/seurat-object/reference/FetchData.html).
  `NULL` (default) uses the same layer selection as
  `FetchData.CellGraph`. On a `CellGraphList`, a missing layer omits
  only that layer's features; metadata and other requested variables are
  still returned. Remaining names are not looked up in `counts` or other
  layers.

- cells:

  Component IDs to fetch. If `NULL`, all loaded
  [`CellGraph`](CellGraph-class.md) objects are used. Unloaded graphs
  raise an error when they are included in `cells`.

- assay:

  Name of assay to fetch layouts from

## Value

A `tbl_df` with columns `x`, `y`, `z` and any requested `vars`.
`add_marker = TRUE` also adds a `marker` column with the marker label of
each node. Methods that extract from multiple components also include a
`component` column.

## Examples

``` r
library(pixelatorR)

se <- ReadPNA_Seurat(minimal_pna_pxl_file(), verbose = FALSE)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpjKKHFf/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
se <- LoadCellGraphs(se, cells = colnames(se)[1], add_layouts = TRUE, verbose = FALSE)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpjKKHFf/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
cg <- CellGraphs(se)[[1]]

# Coordinates only
layout <- FetchLayoutData(cg)

# Include marker counts
layout <- FetchLayoutData(cg, vars = "B2M")

# Include marker labels
layout <- FetchLayoutData(cg, add_marker = TRUE)

# Combine layouts from a CellGraphList
cgl <- CellGraphs(se)
layout <- FetchLayoutData(cgl, cells = colnames(se)[1], vars = "B2M")

# PNAAssay method
layout <- FetchLayoutData(se[["PNA"]], cells = colnames(se)[1], vars = "B2M")

# Seurat method
layout <- FetchLayoutData(se, cells = colnames(se)[1], vars = "B2M")
```
