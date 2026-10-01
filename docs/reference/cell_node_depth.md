# Set the cell plot node depth scale

Tunes depth sizing for ggplot projections. Apparent point diameter
scales with inverse camera distance so that a point twice as far from
the camera is drawn at half the diameter. The point at mean depth keeps
the constant node size from `cell_plot(size = <number>)` (or `1` when no
constant is set).

## Usage

``` r
cell_node_depth(object, focal_distance = 1.5)
```

## Arguments

- object:

  A `cell_plot` recipe.

- focal_distance:

  Positive finite focal distance in the units of the mapped depth
  column. Defaults to `1.5`.

## Value

A modified `cell_plot` recipe.

## Details

Depth sizing is on by default because `depth` is mapped to the `z`
column. It is turned off with `cell_plot(depth = NULL)`, which also
makes this modifier an error. `focal_distance` controls the strength of
the effect. Smaller values exaggerate size differences by depth; larger
values flatten them.

Depth sizing scales one node size, so this modifier cannot be combined
with [`cell_node_scale_size()`](cell_node_scale_size.md). If a size
column is mapped, depth sizing is ignored.
[`cell_plot_interactive()`](cell_plot_interactive.md) and
[`cell_plot_rgl()`](cell_plot_rgl.md) ignore this modifier.

## See also

[`cell_plot()`](cell_plot.md)

Other cell-plot-modifiers: [`cell_annotation()`](cell_annotation.md),
[`cell_grid()`](cell_grid.md),
[`cell_illuminate()`](cell_illuminate.md),
[`cell_node_scale_alpha()`](cell_node_scale_alpha.md),
[`cell_node_scale_color()`](cell_node_scale_color.md),
[`cell_node_scale_size()`](cell_node_scale_size.md),
[`cell_theme()`](cell_theme.md)

## Examples

``` r
se <- ReadPNA_Seurat(minimal_pna_pxl_file())
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpmC3mql/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> ✔ Created a <Seurat> object with 5 cells and 158 targeted surface proteins
se <- LoadCellGraphs(se, cells = colnames(se)[4], verbose = FALSE) |>
  ComputeLayout(layout_method = "spectral")
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpmC3mql/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> ℹ Computing layouts for 1 graphs

cell_graph <- CellGraphs(se)[[4]]

layout_data <- FetchLayoutData(cell_graph, vars = "CD82", layout_method = "spectral_3d") |>
  # Downsample to speed up tests
  dplyr::slice_sample(n = 5000)

cell_plot(layout_data) |>
  cell_node_depth(focal_distance = 5)

```
