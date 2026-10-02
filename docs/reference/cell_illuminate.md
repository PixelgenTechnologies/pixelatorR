# Add illumination to a cell plot

Records geometry-based illumination parameters. During
[`build_cell_plot()`](build_cell_plot.md), directional light, radial
volume shading, and ambient occlusion are calculated independently
within each non-empty panel after color-scale training. The resulting
mask changes rendered point colors, while the unilluminated scale
remains the source for the legend.

## Usage

``` r
cell_illuminate(
  object,
  clamp_quantiles = c(0.01, 0.95),
  directional_light_weight = 0.7,
  volume_shading_weight = 0.5,
  ambient_occlusion_weight = 1,
  ambient_occlusion_k = 20,
  ambient_intensity = 0.3,
  saturation_boost = 0.6,
  shadow_colors = NULL,
  light_direction = c(-3, 2, 3),
  lock_light = FALSE
)
```

## Arguments

- object:

  A `cell_plot` recipe.

- clamp_quantiles:

  Two ordered quantiles between zero and one used to clamp illumination.

- directional_light_weight, volume_shading_weight,
  ambient_occlusion_weight:

  Non-negative weights for light from `light_direction`, radial shading
  from the origin, and local neighbor density, respectively. At least
  one weight must be positive.

- ambient_occlusion_k:

  Positive integer number of nearest neighbors used for ambient
  occlusion.

- ambient_intensity:

  Minimum brightness in fully shadowed regions, between zero and one.

- saturation_boost:

  Non-negative saturation increase in shadowed regions. Ignored when
  `shadow_colors` is supplied.

- shadow_colors:

  Optional vector of fully opaque colors used to tint shadows. When
  `NULL`, illumination is applied in HSV color space.

- light_direction:

  Three finite numeric values giving the x, y, and z direction of the
  light. The vector is normalized before it is stored. The default is a
  late afternoon sun, 45 degrees to the left of the camera and 25
  degrees above the horizon. A light pointing straight at the camera, as
  in `light_direction = c(0, 0, 1)`, instead shades nodes by their
  depth.

- lock_light:

  Whether an animation keeps `light_direction` fixed while rotating the
  cell and recomputes illumination for every frame. This setting is only
  consumed by animation rendering; static and interactive plots use
  illumination calculated from the unrotated coordinates.

## Value

A modified `cell_plot` recipe.

## Details

Illumination masks are normalized within each panel. The requested
number of ambient-occlusion neighbors is capped to the available points.
Singleton panels and panels with a constant mask retain their original
colors. Colors assigned to missing mapped values are not illuminated.
Locked-light animations reuse rotation-invariant ambient occlusion and,
when rotating around the coordinate origin, radial volume shading.

## See also

[`cell_plot_animate()`](cell_plot_animate.md),
[`cell_coord_rotate()`](cell_coord_rotate.md)

Other cell-plot-modifiers: [`cell_annotation()`](cell_annotation.md),
[`cell_grid()`](cell_grid.md),
[`cell_node_depth()`](cell_node_depth.md),
[`cell_node_scale_alpha()`](cell_node_scale_alpha.md),
[`cell_node_scale_color()`](cell_node_scale_color.md),
[`cell_node_scale_size()`](cell_node_scale_size.md),
[`cell_theme()`](cell_theme.md)

## Examples

``` r
se <- ReadPNA_Seurat(minimal_pna_pxl_file())
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpjKKHFf/duckdb
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
#> ℹ /tmp/RtmpjKKHFf/duckdb
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
  cell_illuminate()

```
