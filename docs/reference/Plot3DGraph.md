# Plot 3D graph layouts

**\[deprecated\]**

## Usage

``` r
Plot3DGraph(
  object,
  cell_id,
  marker = NULL,
  assay = NULL,
  layout_method = c("cpmds_3d", "wpmds_3d", "pmds_3d"),
  project = FALSE,
  aspectmode = c("data", "cube"),
  colors = c("lightgrey", "mistyrose", "red", "darkred"),
  showgrid = TRUE,
  log_scale = TRUE,
  node_size = 2,
  show_Bnodes = FALSE,
  ...
)
```

## Arguments

- object:

  A `Seurat` object

- cell_id:

  ID of component to visualize

- marker:

  Name of marker to color the nodes by

- assay:

  Name of assay to pull data from

- layout_method:

  Select appropriate layout previously computed with
  [`ComputeLayout`](ComputeLayout.md)

- project:

  Project the nodes onto a sphere. Default FALSE

- aspectmode:

  Set aspect ratio to one of "data" or "cube". If "cube", this scene's
  axes are drawn as a cube, regardless of the axes' ranges. If "data",
  this scene's axes are drawn in proportion with the axes' ranges.

  Default "data"

- colors:

  Color the nodes expressing a marker. Must be a character vector with
  at least two color names.

- showgrid:

  Show the grid lines. Default TRUE

- log_scale:

  Convert node counts to log-scale with `logp`

- node_size:

  Size of nodes

- show_Bnodes:

  Should B nodes be included in the visualization? This option is only
  applicable to bipartite graphs.

- ...:

  Additional parameters passed to `plot_ly`

## Value

A interactive 3D plot of a component graph layout as a `plotly` object

## Details

Deprecated. Use [`cell_plot()`](cell_plot.md) instead.

Plot a 3D component graph layout computed with
[`ComputeLayout`](ComputeLayout.md) and color nodes by a marker.

## See also

[`cell_plot()`](cell_plot.md)

## Examples

``` r
library(pixelatorR)

# Use cell_plot() instead
if (FALSE) { # \dontrun{
# MPX
pxl_file <- minimal_mpx_pxl_file()
seur <- ReadMPX_Seurat(pxl_file)
seur <- LoadCellGraphs(seur, cells = colnames(seur)[5])
seur <- ComputeLayout(seur, layout_method = "wpmds", dim = 3, pivots = 50)
Plot3DGraph(seur, cell_id = colnames(seur)[5], marker = "CD50", layout_method = "wpmds_3d")

# PNA
pxl_file <- minimal_pna_pxl_file()
seur <- ReadPNA_Seurat(pxl_file)
seur <- LoadCellGraphs(seur, cells = colnames(seur)[1], add_layouts = TRUE)
Plot3DGraph(seur, cell_id = colnames(seur)[1], marker = "CD16", layout_method = "wpmds_3d")
} # }
```
