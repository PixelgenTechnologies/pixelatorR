# Build a cell plot recipe

Resolves mappings, constants, scales, illumination, and display defaults
without drawing a plot. The returned object has the same top-level
fields as the input recipe, with data ordered by the `arrange` mapping.
Rows with missing values (`NA`) in mapped `size` or `alpha` columns are
removed with a warning before color illumination is calculated. Color,
size, and alpha resolve to one value each when the aesthetic is constant
or unmapped. Named categorical size and alpha scales assign one resolved
value per level. Illumination is calculated after color-scale training
and stored separately from the unilluminated resolved colors used by
legends. Other projection-specific transformations, such as apparent
size from depth, are applied by renderers rather than during the build.

## Usage

``` r
build_cell_plot(object)
```

## Arguments

- object:

  A `cell_plot` recipe.

## Value

A `cell_plot_built` object containing resolved plot data and rendering
instructions.

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

built <- cell_plot(layout_data, color = "CD82") |>
  build_cell_plot()
```
