# Compute a graph layout

Computes graph layouts for component graphs.

## Usage

``` r
ComputeLayout(object, ...)

# S3 method for class 'tbl_graph'
ComputeLayout(
  object,
  layout_method = c("pmds", "wpmds"),
  dim = 3,
  normalize_layout = FALSE,
  project_on_unit_sphere = FALSE,
  pivots = 100,
  seed = 123,
  custom_layout_function = NULL,
  custom_layout_function_args = NULL,
  ...
)

# S3 method for class 'CellGraph'
ComputeLayout(
  object,
  layout_method = c("pmds", "wpmds"),
  layout_name = NULL,
  dim = 3,
  normalize_layout = FALSE,
  project_on_unit_sphere = FALSE,
  pivots = 100,
  seed = 123,
  custom_layout_function = NULL,
  custom_layout_function_args = NULL,
  ...
)

# S3 method for class 'MPXAssay'
ComputeLayout(
  object,
  layout_method = c("pmds", "wpmds"),
  layout_name = NULL,
  dim = 3,
  normalize_layout = FALSE,
  project_on_unit_sphere = FALSE,
  pivots = 100,
  seed = 123,
  verbose = TRUE,
  custom_layout_function = NULL,
  custom_layout_function_args = NULL,
  cl = NULL,
  ...
)

# S3 method for class 'CellGraphAssay'
ComputeLayout(
  object,
  layout_method = c("pmds", "wpmds"),
  layout_name = NULL,
  dim = 3,
  normalize_layout = FALSE,
  project_on_unit_sphere = FALSE,
  pivots = 100,
  seed = 123,
  verbose = TRUE,
  custom_layout_function = NULL,
  custom_layout_function_args = NULL,
  cl = NULL,
  ...
)

# S3 method for class 'CellGraphAssay5'
ComputeLayout(
  object,
  layout_method = c("pmds", "wpmds"),
  layout_name = NULL,
  dim = 3,
  normalize_layout = FALSE,
  project_on_unit_sphere = FALSE,
  pivots = 100,
  seed = 123,
  verbose = TRUE,
  custom_layout_function = NULL,
  custom_layout_function_args = NULL,
  cl = NULL,
  ...
)

# S3 method for class 'PNAAssay'
ComputeLayout(
  object,
  layout_method = c("pmds", "wpmds"),
  layout_name = NULL,
  dim = 3,
  normalize_layout = FALSE,
  project_on_unit_sphere = FALSE,
  pivots = 100,
  seed = 123,
  verbose = TRUE,
  custom_layout_function = NULL,
  custom_layout_function_args = NULL,
  cl = NULL,
  ...
)

# S3 method for class 'PNAAssay5'
ComputeLayout(
  object,
  layout_method = c("pmds", "wpmds"),
  layout_name = NULL,
  dim = 3,
  normalize_layout = FALSE,
  project_on_unit_sphere = FALSE,
  pivots = 100,
  seed = 123,
  verbose = TRUE,
  custom_layout_function = NULL,
  custom_layout_function_args = NULL,
  cl = NULL,
  ...
)

# S3 method for class 'Seurat'
ComputeLayout(
  object,
  assay = NULL,
  layout_method = c("pmds", "wpmds"),
  layout_name = NULL,
  dim = 3,
  normalize_layout = FALSE,
  project_on_unit_sphere = FALSE,
  pivots = 100,
  seed = 123,
  verbose = TRUE,
  custom_layout_function = NULL,
  custom_layout_function_args = NULL,
  cl = NULL,
  ...
)
```

## Arguments

- object:

  An object

- ...:

  Additional parameters passed to other methods

- layout_method:

  The method for calculating the graph layout: weighted Pivot MDS
  ("wpmds") or Pivot MDS (pmds)

- dim:

  An integer specifying the dimensions of the layout (2 or 3). Note that
  for "pmds" and "wpmds", the x and y coordinates will be identical
  regardless if dim is set to 2 or 3.

- normalize_layout:

  Logical specifying whether the coordinate system should be centered at
  origo and the coordinates scaled such that their median length
  (euclidean norm) is 1.

- project_on_unit_sphere:

  Should the resulting layout be projected onto a unit sphere?

- pivots:

  Only used for "wpmds" and "pmds" layout algorithms. See
  `?layout_with_pmds` for details

- seed:

  Set seed for reproducibility

- custom_layout_function:

  A custom function for layout computation. The function should take a
  `tbl_graph` object and a `dim` value as input and return a matrix of
  dimensions NxD, where N is the number of nodes and D is equal to
  `dim`. Note that this will override the `layout_method`.

- custom_layout_function_args:

  A list of arguments passed to `custom_layout_function`. The `dim` is
  automatically passed to `custom_layout_function` and should not be
  included in `custom_layout_function_args`.

- layout_name:

  The name of the computed layout. If this name is not given, the
  `layout_method` will be used as the name. If `dim = 3`, a suffix of
  "\_3d" will be added to the layout name. This behavior is ignored if
  `layout_name` is provided.

- verbose:

  Print messages

- cl:

  An integer to indicate number of child-processes (integer values are
  ignored on Windows) for parallel evaluations. See Details on
  performance in the documentation for `pbapply`. The default is NULL,
  which means that no parallelization is used.

- assay:

  Name of assay to compute layouts for

## Value

An object containing a graph layout

## See also

[`center_layout_coordinates()`](layout_coordinates_utils.md) for
centering layout coordinates,
[`normalize_layout_coordinates()`](layout_coordinates_utils.md) for
normalizing layout coordinates and
[`project_layout_coordinates_on_unit_sphere()`](layout_coordinates_utils.md)
for projecting layout coordinates onto a unit sphere.

## Examples

``` r
library(pixelatorR)
library(dplyr)

pxl_file <- minimal_mpx_pxl_file()

# Load example data
seur <- ReadMPX_Seurat(pxl_file)
#> ✔ Created a 'Seurat' object with 5 cells and 80 targeted surface proteins
#> ! Failed to remove temporary file C:/Users/max/AppData/Local/Temp/RtmpmOhqBt/file62e8357c7ed9.h5ad

# Load 1 cellgraph
seur <- LoadCellGraphs(seur,
  cells = colnames(seur)[1],
  force = TRUE
)
#> →    Loading CellGraphs for 1 cells from sample 1
#> ✔ Successfully loaded 1 CellGraph object(s).

# Get CellGraph
cg <- CellGraphs(seur)[[colnames(seur)[1]]]

# Get tbl_graph object
tbl_graph <- slot(cg, name = "cellgraph")

# Compute layout for a tbl_graph
layout <- ComputeLayout(tbl_graph, layout_method = "wpmds")
layout %>% head()
#> # A tibble: 6 × 3
#>            x          y        z
#>        <dbl>      <dbl>    <dbl>
#> 1 -13803157.   3752485. 1268494.
#> 2 -16068734. -11299961. 3206249.
#> 3   2721246.  -4099649. 1741382.
#> 4  13492167.  -7284634. 4530210.
#> 5 -11602483.   7420399. 3358749.
#> 6 -13331214.  -8361592. 4249060.


# Compute layout for a CellGraph
cg <- ComputeLayout(cg, layout_method = "wpmds")


# Compute layout for a CellGraphAssay
cg_assay <- ComputeLayout(seur[["mpxCells"]], layout_method = "wpmds")
#> ℹ Computing layouts for 1 graphs
```
