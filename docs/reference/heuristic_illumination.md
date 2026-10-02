# Compute heuristic illumination for 3D layouts

Combines three simple lighting heuristics for 3D coordinates: (1)
directional light along `light_direction` (default: the positive
z-axis), (2) radial volume shading from the origin, and (3) ambient
occlusion approximated from mean distance to nearest neighbors.

## Usage

``` r
heuristic_illumination(
  layout,
  clamp_quantiles = c(0.01, 0.95),
  directional_light_weight = 0.7,
  volume_shading_weight = 0.5,
  ambient_occlusion_weight = 1,
  ambient_occlusion_k = 20,
  normalize_weights = TRUE,
  light_direction = c(0, 0, 1)
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

- light_direction:

  Numeric vector of length 3 in layout `(x, y, z)` coordinates giving
  the directional light axis. Internally normalized to unit length. The
  directional term is the projection of each point onto this unit
  vector, then rescaled to `[0, 1]`. Default: `c(0, 0, 1)` (positive
  z-axis). The zero vector and non-finite values are rejected.

## Value

A numeric vector of illumination values (length `nrow(layout)`). Higher
values indicate stronger illumination.

## Details

Directional lighting is defined in layout `(x, y, z)` coordinates, not
camera coordinates. Interactive cameras will not re-light a scene unless
a renderer recomputes the illumination mask.

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

# Use cell_plot() and cell_plot_animate() to render a rotating layout
if (FALSE) { # \dontrun{
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
} # }
```
