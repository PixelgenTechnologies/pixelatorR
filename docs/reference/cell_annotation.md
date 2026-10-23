# Annotate a cell plot

Sets title, subtitle, and color legend title text shared by cell plot
renderers. Annotation content is equivalent across backends, although
placement and typography may differ. When `legend_title` is omitted, the
mapped color column name is used.

## Usage

``` r
cell_annotation(object, title = NULL, subtitle = NULL, legend_title = NULL)
```

## Arguments

- object:

  A `cell_plot` recipe.

- title, subtitle, legend_title:

  Optional scalar character strings.

## Value

A modified `cell_plot` recipe.

## See also

Other cell-plot-modifiers: [`cell_grid()`](cell_grid.md),
[`cell_illuminate()`](cell_illuminate.md),
[`cell_node_depth()`](cell_node_depth.md),
[`cell_node_scale_alpha()`](cell_node_scale_alpha.md),
[`cell_node_scale_color()`](cell_node_scale_color.md),
[`cell_node_scale_size()`](cell_node_scale_size.md),
[`cell_theme()`](cell_theme.md)

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

cell_plot(layout_data, color = CD82) |>
  cell_annotation(
    title = "CD82 distribution on a cell",
    subtitle = "Positive nodes colored in red",
    legend_title = "CD82 nodes"
  )

```
