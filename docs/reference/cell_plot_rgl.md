# Render a cell plot with rgl

Builds a [`cell_plot()`](cell_plot.md) recipe and draws it as an
interactive native 3D scatter using
[rgl::rgl](https://dmurdoch.github.io/rgl/dev/reference/rgl-package.html).
Node sizes are converted from backend-neutral relative units to rgl
point diameters in pixels. Continuous sizes are grouped into a bounded
number of pixel-size bins to keep the scene responsive. Panel grids use
[`rgl::layout3d()`](https://dmurdoch.github.io/rgl/dev/reference/mfrow3d.html)
with shared mouse control among data panels. Facet labels are drawn in
dedicated themeable strip regions along the top (columns) and side
(rows), so points cannot cover them. Plot titles and color legends also
use reserved full-window chrome regions. Rotating and scroll-zooming
work with the pointer anywhere in the window, including over the strips,
title, and legend, and always drive every data panel together. Panel
grids support at most 10 rows and 20 columns.

## Usage

``` r
cell_plot_rgl(object)
```

## Arguments

- object:

  A `cell_plot` recipe.

## Value

The rgl device id, returned invisibly.

## Details

Occlusion follows the scene camera, so markers closer to the current
viewpoint appear in front. The rgl backend does not provide hover labels
or interactive legend filtering. The `arrange` and `depth` mappings and
[`cell_coord_rotate()`](cell_coord_rotate.md) are ignored. Markers use
the rendered colors resolved by
[`build_cell_plot()`](build_cell_plot.md), including illumination when
requested, while legends are drawn from the unilluminated color scale
metadata as static 2D overlays. Numeric color mappings use a continuous
colorbar; categorical mappings use a discrete legend.

Unlike [`cell_plot_interactive()`](cell_plot_interactive.md), this
renderer opens an rgl device rather than returning an htmlwidget.
Legends are drawn with
[`rgl::bgplot3d()`](https://dmurdoch.github.io/rgl/dev/reference/bgplot3d.html)
and do not support Plotly-style interactive legend filtering.

[`rgl::bgplot3d()`](https://dmurdoch.github.io/rgl/dev/reference/bgplot3d.html)
renders text into a bitmap sized for the window it was drawn in. Titles,
facet strips, and legends are redrawn automatically at the new size
shortly after the window stops changing.

## See also

[`cell_plot()`](cell_plot.md),
[`cell_plot_interactive()`](cell_plot_interactive.md)

## Examples

``` r
if (FALSE) { # interactive()
se <- ReadPNA_Seurat(minimal_pna_pxl_file())
se <- LoadCellGraphs(se, cells = colnames(se)[4], verbose = FALSE) |>
  ComputeLayout(layout_method = "spectral")

cell_graph <- CellGraphs(se)[[4]]

layout_data <- FetchLayoutData(cell_graph, vars = "CD82", layout_method = "spectral_3d")

cell_plot(layout_data, color = CD82) |>
  cell_plot_rgl()
}
```
