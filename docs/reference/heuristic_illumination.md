# Compute heuristic illumination for 3D layouts

Combines three simple lighting heuristics for 3D coordinates: (1)
directional light from the positive z-axis, (2) radial volume shading
from the origin, and (3) ambient occlusion approximated from mean
distance to nearest neighbors.

## Usage

``` r
heuristic_illumination(
  layout,
  clamp_quantiles = c(0.01, 0.95),
  directional_light_weight = 0.7,
  volume_shading_weight = 0.5,
  ambient_occlusion_weight = 1,
  ambient_occlusion_k = 20,
  normalize_weights = TRUE
)
```

## Arguments

- layout:

  A data frame or tibble with numeric columns `x`, `y`, and `z`.

- clamp_quantiles:

  Numeric vector of length 2 in `[0, 1]`. Illumination is clamped to
  these quantiles to reduce outlier influence. Default: `c(0.01, 0.95)`.

- directional_light_weight:

  Non-negative numeric scalar. Weight for directional light component.
  Default: `0.7`.

- volume_shading_weight:

  Non-negative numeric scalar. Weight for radial volume shading
  component. Default: `0.5`.

- ambient_occlusion_weight:

  Non-negative numeric scalar. Weight for ambient occlusion component.
  Default: `1`.

- ambient_occlusion_k:

  Positive integer. Number of nearest neighbors used for ambient
  occlusion approximation. Default: `20`.

- normalize_weights:

  Logical; if `TRUE`, weights are normalized to sum to 1. Default:
  `TRUE`.

## Value

A numeric vector of illumination values (length `nrow(layout)`). Higher
values indicate stronger illumination.

## Examples

``` r

library(dplyr)
set.seed(1)

# Here we simulate some 3D coordinates with a roughly spherical distribution
n <- 20000
n_surface <- 19000
n_interior <- 1000

# Surface points: normalize to unit sphere, add small Gaussian noise
xyz_surface <- matrix(rnorm(n_surface * 3), ncol = 3)
xyz_surface <- xyz_surface / sqrt(rowSums(xyz_surface^2)) # project to unit sphere
xyz_surface <- xyz_surface + matrix(rnorm(n_surface * 3, sd = 0.05), ncol = 3)

# Interior points: uniform in ball via rejection sampling
xyz_interior <- matrix(rnorm(n_interior * 3), ncol = 3)
radii <- runif(n_interior)^(1 / 3) # cube root for uniform volume distribution
xyz_interior <- xyz_interior / sqrt(rowSums(xyz_interior^2)) * radii * 0.8

layout <- tibble::tibble(
  x = c(xyz_surface[, 1], xyz_interior[, 1]),
  y = c(xyz_surface[, 2], xyz_interior[, 2]),
  z = c(xyz_surface[, 3], xyz_interior[, 3])
)
illum <- heuristic_illumination(layout)

# Create a temporary GIF file and render a rotating layout
# using the computed illumination as node values
temp_gif <- fs::file_temp(ext = ".gif")
render_rotating_layout(
  data = layout %>%
    mutate(node_val = illum),
  pt_size = 0.8,
  width = 740,
  height = 650,
  colors = PixelgenGradient(100, "NaturalBlue"),
  file = temp_gif,
  max_degree = 30,
  frames = 20,
  delay = 1 / 20,
  res = 100,
  boomerang = TRUE,
  show_first_frame = FALSE
)
#> ! The following parameters are not supported with graphics_use = 'base':
#> ggplot_theme, label_grid_axes, title
#> Titles, text annotations and color bar will be missing from the plot.
#> 
#> ── Rendering frames... 
#> 
#> ── Creating video... 
#> ✔ The video has been saved to:
#> /var/folders/gw/bdcqhnvs0m9gs_mq8n51jtbc0000gn/T/RtmpemhbYW/filea21a4ef05ae5.gif
```
