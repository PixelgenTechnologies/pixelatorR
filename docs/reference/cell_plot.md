# Create a cell plot recipe

Creates a cell-layout visualization recipe. The recipe stores data,
column mappings, and optional rendering instructions without drawing a
plot. Add instructions with the `cell_*()` modifier functions. Printing
the recipe builds and renders a ggplot.
[`cell_plot_interactive()`](cell_plot_interactive.md) draws the same
recipe as a 3D Plotly widget. [`cell_plot_rgl()`](cell_plot_rgl.md)
draws it as a native rgl scene.
[`cell_plot_animate()`](cell_plot_animate.md) encodes a rotating GIF or
video after [`cell_coord_rotate()`](cell_coord_rotate.md).
[`summary()`](https://rdrr.io/r/base/summary.html) inspects the recipe
without drawing.

## Usage

``` r
cell_plot(
  data,
  x = x,
  y = y,
  z = z,
  color = NULL,
  size = NULL,
  alpha = NULL,
  arrange = z,
  depth = z,
  illumination_mask = NULL
)
```

## Arguments

- data:

  A non-empty tibble containing cell-layout data.

- x, y:

  Column mappings for the horizontal and vertical coordinates. Bare
  column names and character column names are supported.

- z:

  Column mapping for the third coordinate. Defaults to the `z` column.

- color:

  Optional column mapping or constant color for every node. Character,
  factor, numeric, and integer columns are supported, where character
  and factor columns are treated as categorical and numeric and integer
  columns as continuous. A character value is read as a column name when
  `data` has that column and otherwise as one fully opaque color.

- size, alpha:

  Optional column mappings or constant values for node size and alpha.
  Character, factor, numeric, and integer columns are supported, where
  character and factor columns are treated as categorical levels. One
  finite number sets a constant relative size or a constant alpha
  between zero and one. Category-specific sizes and alphas are set with
  [`cell_node_scale_size()`](cell_node_scale_size.md) and
  [`cell_node_scale_alpha()`](cell_node_scale_alpha.md).

- arrange:

  An optional column mapping used to order points in projected plots.
  [`cell_plot_interactive()`](cell_plot_interactive.md) and
  [`cell_plot_rgl()`](cell_plot_rgl.md) ignore this mapping because
  occlusion follows the scene camera. Defaults to the mapped `z` column.

- depth:

  An optional column mapping used for ggplot depth sizing. Defaults to
  the mapped `z` column, so depth sizing is on unless `depth = NULL`
  turns it off. Must not be the same column as `x`, `y`, or `size`.
  Plotly and rgl ignore this mapping.

- illumination_mask:

  Optional logical column mapping. Illumination is applied only to rows
  where this column is `TRUE`. Rows where it is `FALSE` retain their
  unilluminated colors. Missing values are not allowed.

## Value

A `cell_plot` recipe.

## Details

A `cell_plot` is an S3 list with a fixed set of fields:

- `data`: the visualization tibble

- `mapping`: named column mappings (`x`, `y`, `z`, `depth`, `color`,
  `size`, `alpha`, `illumination_mask`, `arrange`)

- `constant`: single `color`, `size`, and `alpha` values shared by every
  point, `NULL` when the aesthetic is mapped or left at its default

- `grid`, `color`, `size`, `depth`, `alpha`, `theme`, `annotation`,
  `illuminate`, `coord`: instruction slots, `NULL` until a modifier
  writes them

Coordinate defaults are `x`, `y`, and `z`. `arrange` and `depth` default
to the mapped `z` column. Pass `depth = NULL` to turn off depth sizing.
Color, size, and alpha default to no mapping.

`color`, `size`, and `alpha` accept either a column mapping or one
constant value shared by every point, such as `color = "red"`,
`size = 3`, or `alpha = 0.5`. A constant aesthetic has no scale, so the
matching `cell_node_scale_*()` modifier cannot be used with it. Without
a mapping and without a constant, nodes are `gray90` with relative size
`1` and alpha `1`.

When a size column is mapped, the default output range is `c(2, 6)`
unless [`cell_node_scale_size()`](cell_node_scale_size.md) overrides it.
For a categorical size or alpha mapping,
[`cell_node_scale_size()`](cell_node_scale_size.md) and
[`cell_node_scale_alpha()`](cell_node_scale_alpha.md) also accept a
named vector of per-level values, such as
`sizes = c(a = 1, b = 2, c = 3)`. Color scale behavior is controlled by
[`cell_node_scale_color()`](cell_node_scale_color.md). Color scales are
trained once across all panels, and continuous scales span the observed
range unless `limits` are supplied.

## See also

[`cell_node_scale_size()`](cell_node_scale_size.md),
[`cell_node_depth()`](cell_node_depth.md),
[`cell_coord_rotate()`](cell_coord_rotate.md),
[`cell_plot_animate()`](cell_plot_animate.md),
[`cell_plot_interactive()`](cell_plot_interactive.md),
[`cell_plot_rgl()`](cell_plot_rgl.md)

## Examples

``` r
# Plot a spectral layout of a cell from the example data
se <- ReadPNA_Seurat(minimal_pna_pxl_file())
#> ✔ Created a <Seurat> object with 5 cells and 158 targeted surface proteins
se <- LoadCellGraphs(se, cells = colnames(se)[4], verbose = FALSE) |>
  ComputeLayout(layout_method = "spectral")
#> ℹ Computing layouts for 1 graphs

cell_graph <- CellGraphs(se)[[4]]

layout_data <- FetchLayoutData(cell_graph, vars = "CD82", layout_method = "spectral_3d")

cell_plot(layout_data, color = CD82) |>
  cell_node_scale_color(colors = c("lightgrey", "red")) |>
  cell_annotation(title = "Spectral layout", subtitle = colnames(se)[1])

```
