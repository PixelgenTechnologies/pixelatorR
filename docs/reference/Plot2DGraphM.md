# Plot multiple markers on multiple graphs

In contrast to [`Plot2DGraph`](Plot2DGraph.md), which only draw 1 marker
at the time, this function makes it possible to arrange plots into a
grid with markers in rows and components in columns. The color scales
are fixed for each marker so that their limits are the same across all
components.

## Usage

``` r
Plot2DGraphM(
  object,
  cells,
  markers,
  assay = NULL,
  layout_method = c("wpmds_3d", "pmds_3d", "wpmds", "pmds"),
  colors = c("lightgrey", "mistyrose", "red", "darkred"),
  map_nodes = TRUE,
  map_edges = FALSE,
  log_scale = TRUE,
  node_size = 0.5,
  edge_width = 0.3,
  show_Bnodes = TRUE,
  titles = NULL,
  titles_theme = NULL,
  titles_size = 10,
  titles_col = "black",
  ...
)
```

## Arguments

- object:

  A `Seurat` object

- cells:

  A character vector with cell IDs

- markers:

  A character vector with marker names

- assay:

  Name of assay to pull data from

- layout_method:

  Select appropriate layout previously computed with
  [`ComputeLayout`](ComputeLayout.md)

- colors:

  A character vector of colors to color marker counts by

- map_nodes, map_edges:

  Should nodes and/or edges be mapped? Note that component graphs can
  have \>100k edges which can be very slow to draw.

- log_scale:

  Convert node counts to log-scale with `log1p`. This parameter is
  ignored for PNA graphs.

- node_size:

  Size of nodes

- edge_width:

  Set the width of the edges if `map_edges = TRUE`

- show_Bnodes:

  Should B nodes be included in the visualization? This option is only
  applicable to bipartite MPX graphs. Note that by removing the B nodes,
  all edges are removed from the graph and hence, `map_edges` will have
  no effect.

- titles:

  A named character vector with optional titles. The names of `titles`
  should match `cells`

- titles_theme:

  A `theme` used to style the titles

- titles_size:

  The size of the text in the plot titles

- titles_col:

  The color of the plot titles

- ...:

  Parameters passed to Plot2DGraph

## Value

A `patchwork` object

## See also

[`Plot2DGraph()`](Plot2DGraph.md)

## Examples

``` r
library(pixelatorR)

# MPX
pxl_file <- minimal_mpx_pxl_file()
seur <- ReadMPX_Seurat(pxl_file)
#> ! Failed to remove temporary dir C:/Users/max/AppData/Local/Temp/RtmpInAYcT/dir7de03d421221
#> ! Failed to remove temporary dir C:/Users/max/AppData/Local/Temp/RtmpInAYcT/dir7de03242294
#> ! Failed to remove temporary file C:/Users/max/AppData/Local/Temp/RtmpInAYcT/file7de07e4bf38.h5ad
seur <- LoadCellGraphs(seur, load_as = "Anode")
seur <- ComputeLayout(seur, layout_method = "pmds", dim = 2)
Plot2DGraphM(seur, cells = colnames(seur)[2:3], layout_method = "pmds", markers = c("CD20", "CD4"))
#> ✖ 'CD20' is missing from node count matrix for component RCVCMP0000487


# PNA
pxl_file <- minimal_pna_pxl_file()
seur <- ReadPNA_Seurat(pxl_file)
seur <- LoadCellGraphs(seur, cells = colnames(seur)[2:3], add_layouts = TRUE)
Plot2DGraphM(seur, cells = colnames(seur)[2:3], markers = c("CD20", "CD4"))

```
