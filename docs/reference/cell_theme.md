# Set the cell plot theme

Sets the backend-neutral appearance shared by cell plot renderers. The
contract contains `background_color`, `strip_background_color`,
`text_color`, and `text_size`. Backends translate these values to their
own theme systems. Backend-specific theme objects are not stored here.

## Usage

``` r
cell_theme(
  object,
  background_color,
  text_color,
  text_size,
  strip_background_color
)
```

## Arguments

- object:

  A `cell_plot` recipe.

- background_color:

  Background color. Defaults to white.

- text_color:

  Text color. Defaults to black.

- text_size:

  Positive text size. Defaults to 11.

- strip_background_color:

  Facet strip background color. Defaults to the ggplot2-like gray
  `"#D9D9D9"`.

## Value

A modified `cell_plot` recipe.

## See also

Other cell-plot-modifiers: [`cell_annotation()`](cell_annotation.md),
[`cell_grid()`](cell_grid.md),
[`cell_illuminate()`](cell_illuminate.md),
[`cell_node_depth()`](cell_node_depth.md),
[`cell_node_scale_alpha()`](cell_node_scale_alpha.md),
[`cell_node_scale_color()`](cell_node_scale_color.md),
[`cell_node_scale_size()`](cell_node_scale_size.md)

## Examples

``` r
se <- ReadPNA_Seurat(minimal_pna_pxl_file())
#> ✔ Created a <Seurat> object with 5 cells and 158 targeted surface proteins
se <- LoadCellGraphs(se, cells = colnames(se)[4], verbose = FALSE) |>
  ComputeLayout(layout_method = "spectral")
#> ℹ Computing layouts for 1 graphs

cell_graph <- CellGraphs(se)[[4]]

layout_data <- FetchLayoutData(cell_graph, vars = "CD82", layout_method = "spectral_3d") |>
  # Downsample to speed up tests
  dplyr::slice_sample(n = 5000)

cell_plot(layout_data) |>
  cell_theme(
    background_color = "navy",
    strip_background_color = "grey30",
    text_color = "white",
    text_size = 12
  )

```
