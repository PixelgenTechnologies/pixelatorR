# Set the cell plot node alpha scale

Defines the output alphas for a column mapped with
`cell_plot(alpha = ...)`. Numeric and integer columns are scaled
continuously. Character and factor columns are treated as categorical
levels.

## Usage

``` r
cell_node_scale_alpha(object, alphas = c(0.2, 1), limits = NULL)
```

## Arguments

- object:

  A `cell_plot` recipe.

- alphas:

  Values between zero and one. Supply two ordered values for a
  continuous or interpolated categorical range, or a named vector of
  per-level values when the mapped column is categorical.

- limits:

  Optional scale limits for a ranged mapping. Supply two ordered numeric
  values. Named categorical alphas cannot be combined with `limits`.

## Value

A modified `cell_plot` recipe.

## Details

A length-2 `alphas` vector is an output range.
[`build_cell_plot()`](build_cell_plot.md) normalizes continuous values,
and interpolated categorical levels, to that interval. Rows with missing
values in a mapped alpha column are dropped at build time with a
warning. A named `alphas` vector assigns one opacity to each category,
for example `alphas = c(a = 0.2, b = 0.5, c = 1)`. Named values must
cover every observed level and are not interpolated. `limits` belong to
the ranged scale and cannot be combined with named alphas.

[`build_cell_plot()`](build_cell_plot.md) resolves alpha between zero
and one. All renderers consume those resolved values directly.

This modifier requires an alpha mapping. One alpha shared by every point
is set with `cell_plot(alpha = <number>)` instead.

## See also

Other cell-plot-modifiers: [`cell_annotation()`](cell_annotation.md),
[`cell_grid()`](cell_grid.md),
[`cell_illuminate()`](cell_illuminate.md),
[`cell_node_depth()`](cell_node_depth.md),
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

cell_plot(layout_data, alpha = CD82) |>
  cell_node_scale_alpha(alphas = c(0.2, 1))


# One alpha per category
layout_data |>
  dplyr::mutate(group = dplyr::if_else(CD82 > median(CD82), "high", "low")) |>
  cell_plot(alpha = group) |>
  cell_node_scale_alpha(alphas = c(low = 0.3, high = 1))


# One alpha for every point
cell_plot(layout_data, alpha = 0.5)

```
