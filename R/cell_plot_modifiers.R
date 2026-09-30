#' Arrange a cell plot in a panel grid
#'
#' `r lifecycle::badge("experimental")`
#'
#' Defines panel rows and columns using existing categorical columns. Named
#' arguments are preferred because they make grid orientation explicit. Facet
#' values may be character, factor, or integer columns. This modifier records
#' mappings only and does not reshape the plot data.
#'
#' @param object A `cell_plot` recipe.
#' @param rows,cols Optional bare or character column names used for panel rows
#' and columns. At least one must be supplied.
#'
#' @return A modified `cell_plot` recipe.
#' @family cell-plot-modifiers
#'
#' @examples
#' library(dplyr)
#' se <- ReadPNA_Seurat(minimal_pna_pxl_file())
#' se <- LoadCellGraphs(se, cells = colnames(se)[3:4], verbose = FALSE) |>
#'   ComputeLayout(layout_method = "cpmds")
#' layout_data <- FetchLayoutData(se, vars = c("CD81", "CD82"), layout_method = "cpmds_3d") |>
#'   # Downsample to speed up tests
#'   group_by(component) |>
#'   slice_sample(n = 5000) |>
#'   tidyr::pivot_longer(cols = c("CD81", "CD82"), names_to = "marker", values_to = "value")
#'
#'
#' cell_plot(layout_data) |>
#'   cell_grid(rows = marker, cols = component)
#'
#' @export
cell_grid <- function(object, rows = NULL, cols = NULL) {
  rows <- .cell_plot_column(rlang::enquo(rows), "rows")
  cols <- .cell_plot_column(rlang::enquo(cols), "cols")
  if (is.null(rows) && is.null(cols)) {
    cli::cli_abort(
      c("i" = "At least one of {.arg rows} or {.arg cols} must be supplied.")
    )
  }

  facet_columns <- unlist(list(rows = rows, cols = cols), use.names = FALSE)
  for (column in facet_columns) {
    pixelatorR:::assert_col_in_data(
      column,
      data = object$data,
      arg_data = "data"
    )
    pixelatorR:::assert_class(
      object$data[[column]],
      classes = c("character", "factor", "integer"),
      arg = column
    )
  }

  specification <- list(rows = rows, cols = cols)
  return(.replace_cell_plot_field(object, "grid", specification))
}

#' Set the cell plot node color scale
#'
#' `r lifecycle::badge("experimental")`
#'
#' Controls how a column mapped with `cell_plot(color = ...)` is represented.
#' Continuous versus categorical scales can be selected explicitly, or inferred
#' from the mapped column. Palette, limits, and missing-value color are recorded
#' here. Color scales are trained once across all panels, and continuous scales
#' span the observed range unless `limits` are supplied. A continuous scale is
#' therefore centered on zero only when `limits` are symmetric around zero, as
#' in `limits = c(-2, 2)`. Palette and missing-value colors must be fully
#' opaque. Node opacity is controlled with [cell_node_scale_alpha()].
#'
#' This modifier requires a color mapping. One color shared by every point is
#' set with `cell_plot(color = "red")` instead.
#'
#' @param object A `cell_plot` recipe.
#' @param colors Optional vector of valid, fully opaque colors used for the
#' scale. When `NULL`, the builder uses its default palette.
#' @param limits Optional scale limits. Supply two ordered numeric values for a
#' continuous scale or a character vector of levels for a categorical scale.
#' @param na_color Fully opaque color used for missing values.
#' @param type Scale type. `"auto"` determines the type from the mapped column.
#'
#' @return A modified `cell_plot` recipe.
#' @family cell-plot-modifiers
#'
#' @examples
#' se <- ReadPNA_Seurat(minimal_pna_pxl_file())
#' se <- LoadCellGraphs(se, cells = colnames(se)[4], verbose = FALSE) |>
#'   ComputeLayout(layout_method = "spectral")
#'
#' cell_graph <- CellGraphs(se)[[4]]
#'
#' layout_data <- FetchLayoutData(cell_graph, vars = "CD82", layout_method = "spectral_3d") |>
#'   # Downsample to speed up tests
#'   dplyr::slice_sample(n = 5000)
#'
#' cell_plot(layout_data, color = CD82) |>
#'   cell_node_scale_color(colors = c("blue", "red"))
#'
#' # One color for every point
#' cell_plot(layout_data, color = "red")
#'
#' @seealso [cell_plot()]
#'
#' @export
cell_node_scale_color <- function(
  object,
  colors = NULL,
  limits = NULL,
  na_color = "grey50",
  type = c("auto", "continuous", "categorical")
) {
  .assert_cell_scale_mapping(
    object,
    aesthetic = "color",
    modifier = "cell_node_scale_color"
  )

  type <- match.arg(type)
  if (!is.null(colors)) {
    pixelatorR:::assert_valid_color(colors, arg = "colors")
    .assert_opaque_colors(colors, arg = "colors")
  }
  pixelatorR:::assert_valid_color(na_color, arg = "na_color")
  pixelatorR:::assert_length(na_color, n = 1, arg_x = "na_color")
  .assert_opaque_colors(na_color, arg = "na_color")
  pixelatorR:::assert_class(
    limits,
    classes = c("numeric", "integer", "character"),
    allow_null = TRUE,
    arg = "limits"
  )

  if (is.numeric(limits)) {
    .assert_ordered_numeric_pair(limits, arg = "limits", limits = c(-Inf, Inf))
  }
  if (is.character(limits)) {
    pixelatorR:::assert_vector(
      limits,
      type = "character",
      n = 1,
      arg = "limits"
    )
    pixelatorR:::assert_unique(limits, arg = "limits")
  }

  color_values <- object$data[[object$mapping$color]]
  scale_type <- .cell_color_scale_type(color_values, type)
  if (scale_type == "continuous" && is.character(limits)) {
    cli::cli_abort(
      c("x" = "Continuous color mappings require numeric {.arg limits}.")
    )
  }
  if (scale_type == "categorical" && is.numeric(limits)) {
    cli::cli_abort(
      c("x" = "Categorical color mappings require character {.arg limits}.")
    )
  }
  if (scale_type == "categorical" && !is.null(names(colors))) {
    color_levels <- limits %||% .cell_categorical_levels(color_values)
    pixelatorR:::assert_x_in_y(
      color_levels,
      names(colors),
      arg_x = "color levels",
      arg_y = "names(colors)"
    )
  }

  specification <- list(
    colors = colors,
    limits = limits,
    na_color = na_color,
    type = type
  )
  return(.replace_cell_plot_field(object, "color", specification))
}

#' Set the cell plot node size scale
#'
#' `r lifecycle::badge("experimental")`
#'
#' Sets the output sizes for a column mapped with `cell_plot(size = ...)`.
#' Numeric and integer columns are scaled continuously. Character and factor
#' columns are treated as categorical levels.
#'
#' A length-2 `sizes` vector is an output range. [build_cell_plot()] normalizes
#' continuous values, and interpolated categorical levels, to that interval.
#' Rows with missing values in a mapped size column are dropped at build time
#' with a warning.
#' A named `sizes` vector assigns one relative size to each category, for
#' example `sizes = c(a = 1, b = 2, c = 3)`. Named values must cover every
#' observed level and are not interpolated. `limits` belong to the ranged
#' scale and cannot be combined with named sizes.
#'
#' Sizes are backend-neutral relative units. Each renderer converts the values
#' to native units: ggplot2 uses millimetres, Plotly and rgl use marker
#' diameters in pixels, and base R uses `cex`. The contract guarantees
#' consistent ordering and relative range, not pixel-perfect physical equality.
#'
#' This modifier requires a size mapping and is optional; the default output
#' range is `c(2, 6)`. One node size shared by every point is set with
#' `cell_plot(size = <number>)` instead. Because [cell_node_depth()] scales one
#' node size by depth, it cannot be combined with this modifier.
#'
#' For projected 3D data, ggplot2 applies perspective sizing from the `depth`
#' mapping only when no size column is mapped. A size mapping always wins.
#' Plotly and rgl ignore `depth`.
#'
#' @param object A `cell_plot` recipe.
#' @param sizes Non-negative relative sizes. Supply two ordered values for a
#' continuous or interpolated categorical range, or a named vector of
#' per-level values when the mapped column is categorical.
#' @param limits Optional scale limits for a ranged mapping. Supply two
#' ordered numeric values. Named categorical sizes cannot be combined with
#' `limits`.
#'
#' @return A modified `cell_plot` recipe.
#' @family cell-plot-modifiers
#'
#' @examples
#' se <- ReadPNA_Seurat(minimal_pna_pxl_file())
#' se <- LoadCellGraphs(se, cells = colnames(se)[4], verbose = FALSE) |>
#'   ComputeLayout(layout_method = "spectral")
#'
#' cell_graph <- CellGraphs(se)[[4]]
#'
#' layout_data <- FetchLayoutData(cell_graph, vars = "CD82", layout_method = "spectral_3d") |>
#'   # Downsample to speed up tests
#'   dplyr::slice_sample(n = 5000)
#'
#' cell_plot(layout_data, size = CD82) |>
#'   cell_node_scale_size(sizes = c(0.1, 3))
#'
#' # One size per category
#' layout_data |>
#'   dplyr::mutate(group = dplyr::if_else(CD82 > median(CD82), "high", "low")) |>
#'   cell_plot(size = group) |>
#'   cell_node_scale_size(sizes = c(low = 1, high = 4))
#'
#' @seealso [cell_plot()], [cell_node_depth()]
#'
#' @export
cell_node_scale_size <- function(object, sizes = c(2, 6), limits = NULL) {
  .assert_cell_numeric_scale_output(
    object,
    aesthetic = "size",
    modifier = "cell_node_scale_size",
    output = sizes,
    output_arg = "sizes",
    output_limits = c(0, Inf),
    limits = limits
  )

  specification <- list(sizes = sizes, limits = limits)
  return(.replace_cell_plot_field(object, "size", specification))
}

#' Set the cell plot node depth scale
#'
#' `r lifecycle::badge("experimental")`
#'
#' Tunes depth sizing for ggplot projections. Apparent point diameter scales
#' with inverse camera distance so that a point twice as far from the camera is
#' drawn at half the diameter. The point at mean depth keeps the constant node
#' size from `cell_plot(size = <number>)` (or `1` when no constant is set).
#'
#' Depth sizing is on by default because `depth` is mapped to the `z` column.
#' It is turned off with `cell_plot(depth = NULL)`, which also makes this
#' modifier an error. `focal_distance` controls the strength of the effect.
#' Smaller values exaggerate size differences by depth; larger values flatten
#' them.
#'
#' Depth sizing scales one node size, so this modifier cannot be combined with
#' [cell_node_scale_size()]. If a size column is mapped, depth sizing is
#' ignored. [cell_plot_interactive()] and [cell_plot_rgl()] ignore this
#' modifier.
#'
#' @param object A `cell_plot` recipe.
#' @param focal_distance Positive finite focal distance in the units of the
#' mapped depth column. Defaults to `1.5`.
#'
#' @return A modified `cell_plot` recipe.
#' @family cell-plot-modifiers
#'
#' @seealso [cell_plot()]
#'
#' @examples
#' se <- ReadPNA_Seurat(minimal_pna_pxl_file())
#' se <- LoadCellGraphs(se, cells = colnames(se)[4], verbose = FALSE) |>
#'   ComputeLayout(layout_method = "spectral")
#'
#' cell_graph <- CellGraphs(se)[[4]]
#'
#' layout_data <- FetchLayoutData(cell_graph, vars = "CD82", layout_method = "spectral_3d") |>
#'   # Downsample to speed up tests
#'   dplyr::slice_sample(n = 5000)
#'
#' cell_plot(layout_data) |>
#'   cell_node_depth(focal_distance = 5)
#'
#' @export
cell_node_depth <- function(object, focal_distance = 1.5) {
  pixelatorR:::assert_single_value(
    focal_distance,
    type = "numeric",
    arg = "focal_distance"
  )
  pixelatorR:::assert_within_limits(
    focal_distance,
    limits = c(0, Inf),
    arg = "focal_distance"
  )
  if (!is.finite(focal_distance) || focal_distance == 0) {
    cli::cli_abort(c("i" = "{.arg focal_distance} must be positive and finite."))
  }

  specification <- list(focal_distance = focal_distance)
  return(.replace_cell_plot_field(object, "depth", specification))
}

.cell_plot_default_focal_distance <- formals(cell_node_depth)$focal_distance

#' Set the cell plot node alpha scale
#'
#' `r lifecycle::badge("experimental")`
#'
#' Defines the output alphas for a column mapped with `cell_plot(alpha = ...)`.
#' Numeric and integer columns are scaled continuously. Character and factor
#' columns are treated as categorical levels.
#'
#' A length-2 `alphas` vector is an output range. [build_cell_plot()]
#' normalizes continuous values, and interpolated categorical levels, to that
#' interval. Rows with missing values in a mapped alpha column are dropped at
#' build time with a warning. A named `alphas` vector assigns one opacity to each category, for
#' example `alphas = c(a = 0.2, b = 0.5, c = 1)`. Named values must cover every
#' observed level and are not interpolated. `limits` belong to the ranged
#' scale and cannot be combined with named alphas.
#'
#' [build_cell_plot()] resolves alpha between zero and one. All renderers
#' consume those resolved values directly.
#'
#' This modifier requires an alpha mapping. One alpha shared by every point is
#' set with `cell_plot(alpha = <number>)` instead.
#'
#' @param object A `cell_plot` recipe.
#' @param alphas Values between zero and one. Supply two ordered values for a
#' continuous or interpolated categorical range, or a named vector of
#' per-level values when the mapped column is categorical.
#' @param limits Optional scale limits for a ranged mapping. Supply two
#' ordered numeric values. Named categorical alphas cannot be combined with
#' `limits`.
#'
#' @return A modified `cell_plot` recipe.
#' @family cell-plot-modifiers
#'
#' @examples
#' se <- ReadPNA_Seurat(minimal_pna_pxl_file())
#' se <- LoadCellGraphs(se, cells = colnames(se)[4], verbose = FALSE) |>
#'   ComputeLayout(layout_method = "spectral")
#'
#' cell_graph <- CellGraphs(se)[[4]]
#'
#' layout_data <- FetchLayoutData(cell_graph, vars = "CD82", layout_method = "spectral_3d") |>
#'   # Downsample to speed up tests
#'   dplyr::slice_sample(n = 5000)
#'
#' cell_plot(layout_data, alpha = CD82) |>
#'   cell_node_scale_alpha(alphas = c(0.2, 1))
#'
#' # One alpha per category
#' layout_data |>
#'   dplyr::mutate(group = dplyr::if_else(CD82 > median(CD82), "high", "low")) |>
#'   cell_plot(alpha = group) |>
#'   cell_node_scale_alpha(alphas = c(low = 0.3, high = 1))
#'
#' # One alpha for every point
#' cell_plot(layout_data, alpha = 0.5)
#'
#' @export
cell_node_scale_alpha <- function(object, alphas = c(0.2, 1), limits = NULL) {
  .assert_cell_numeric_scale_output(
    object,
    aesthetic = "alpha",
    modifier = "cell_node_scale_alpha",
    output = alphas,
    output_arg = "alphas",
    output_limits = c(0, 1),
    limits = limits
  )

  specification <- list(alphas = alphas, limits = limits)
  return(.replace_cell_plot_field(object, "alpha", specification))
}

#' Normalize a direction in cell plot coordinates
#'
#' Validates a three-dimensional direction and scales it to unit length.
#'
#' @param direction A numeric vector with x, y, and z components.
#' @param argument The argument name used in validation errors.
#' @param call The calling environment.
#'
#' @return A numeric vector of length three and unit length.
#'
#' @noRd
.normalize_cell_direction <- function(
  direction,
  argument,
  call = rlang::caller_env()
) {
  pixelatorR:::assert_vector(
    direction,
    type = "numeric",
    n = 3,
    arg = argument,
    call = call
  )
  pixelatorR:::assert_length(
    direction,
    n = 3,
    arg_x = argument,
    call = call
  )
  if (!all(is.finite(direction))) {
    cli::cli_abort(
      c("i" = "{.arg {argument}} must contain only finite values."),
      call = call
    )
  }

  scale <- max(abs(direction))
  if (scale == 0) {
    cli::cli_abort(
      c("i" = "{.arg {argument}} must point in a nonzero direction."),
      call = call
    )
  }

  direction <- direction / scale
  return(direction / sqrt(sum(direction^2)))
}

#' Add illumination to a cell plot
#'
#' `r lifecycle::badge("experimental")`
#'
#' Records geometry-based illumination parameters. During [build_cell_plot()],
#' directional light, radial volume shading, and ambient occlusion are
#' calculated independently within each non-empty panel after color-scale
#' training. The resulting mask changes rendered point colors, while the
#' unilluminated scale remains the source for the legend.
#'
#' Illumination masks are normalized within each panel. The requested number of
#' ambient-occlusion neighbors is capped to the available points. Singleton
#' panels and panels with a constant mask retain their original colors. Colors
#' assigned to missing mapped values are not illuminated. Locked-light
#' animations reuse rotation-invariant ambient occlusion and, when rotating
#' around the coordinate origin, radial volume shading.
#'
#' @param object A `cell_plot` recipe.
#' @param clamp_quantiles Two ordered quantiles between zero and one used to
#' clamp illumination.
#' @param directional_light_weight,volume_shading_weight,ambient_occlusion_weight
#' Non-negative weights for light from `light_direction`, radial shading from
#' the origin, and local neighbor density, respectively. At least one weight
#' must be positive.
#' @param ambient_occlusion_k Positive integer number of nearest neighbors used
#' for ambient occlusion.
#' @param ambient_intensity Minimum brightness in fully shadowed regions,
#' between zero and one.
#' @param saturation_boost Non-negative saturation increase in shadowed regions.
#' Ignored when `shadow_colors` is supplied.
#' @param shadow_colors Optional vector of fully opaque colors used to tint
#' shadows. When `NULL`, illumination is applied in HSV color space.
#' @param light_direction Three finite numeric values giving the x, y, and z
#' direction of the light. The vector is normalized before it is stored. The
#' default is a late afternoon sun, 45 degrees to the left of the camera and 25
#' degrees above the horizon. A light pointing straight at the camera, as in
#' `light_direction = c(0, 0, 1)`, instead shades nodes by their depth.
#' @param lock_light Whether an animation keeps `light_direction` fixed while
#' rotating the cell and recomputes illumination for every frame. This setting
#' is only consumed by animation rendering; static and interactive plots use
#' illumination calculated from the unrotated coordinates.
#'
#' @return A modified `cell_plot` recipe.
#' @family cell-plot-modifiers
#'
#' @seealso [cell_plot_animate()], [cell_coord_rotate()]
#'
#' @examples
#' se <- ReadPNA_Seurat(minimal_pna_pxl_file())
#' se <- LoadCellGraphs(se, cells = colnames(se)[4], verbose = FALSE) |>
#'   ComputeLayout(layout_method = "spectral")
#'
#' cell_graph <- CellGraphs(se)[[4]]
#'
#' layout_data <- FetchLayoutData(cell_graph, vars = "CD82", layout_method = "spectral_3d") |>
#'   # Downsample to speed up tests
#'   dplyr::slice_sample(n = 5000)
#'
#' cell_plot(layout_data) |>
#'   cell_illuminate()
#'
#' @export
cell_illuminate <- function(
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
) {
  pixelatorR:::assert_single_value(
    ambient_intensity,
    type = "numeric",
    arg = "ambient_intensity"
  )
  pixelatorR:::assert_within_limits(
    ambient_intensity,
    limits = c(0, 1),
    arg = "ambient_intensity"
  )
  pixelatorR:::assert_single_value(
    saturation_boost,
    type = "numeric",
    arg = "saturation_boost"
  )
  pixelatorR:::assert_within_limits(
    saturation_boost,
    limits = c(0, Inf),
    arg = "saturation_boost"
  )
  if (!is.finite(saturation_boost)) {
    cli::cli_abort(
      c("i" = "{.arg saturation_boost} must be finite.")
    )
  }
  if (!is.null(shadow_colors)) {
    pixelatorR:::assert_valid_color(shadow_colors, arg = "shadow_colors")
    .assert_opaque_colors(shadow_colors, arg = "shadow_colors")
  }
  light_direction <- .normalize_cell_direction(
    light_direction,
    argument = "light_direction"
  )
  pixelatorR:::assert_single_value(
    lock_light,
    type = "bool",
    arg = "lock_light"
  )

  specification <- list(
    clamp_quantiles = clamp_quantiles,
    directional_light_weight = directional_light_weight,
    volume_shading_weight = volume_shading_weight,
    ambient_occlusion_weight = ambient_occlusion_weight,
    ambient_occlusion_k = ambient_occlusion_k,
    ambient_intensity = ambient_intensity,
    saturation_boost = saturation_boost,
    shadow_colors = shadow_colors,
    light_direction = light_direction,
    lock_light = lock_light
  )
  return(.replace_cell_plot_field(object, "illuminate", specification))
}

#' Set the cell plot theme
#'
#' `r lifecycle::badge("experimental")`
#'
#' Sets the backend-neutral appearance shared by cell plot renderers. The
#' contract contains `background_color`, `strip_background_color`,
#' `text_color`, and `text_size`. Backends translate these values to their own
#' theme systems. Backend-specific theme objects are not stored here.
#'
#' @param object A `cell_plot` recipe.
#' @param background_color Background color. Defaults to white.
#' @param text_color Text color. Defaults to black.
#' @param text_size Positive text size. Defaults to 11.
#' @param strip_background_color Facet strip background color. Defaults to the
#'   ggplot2-like gray `"#D9D9D9"`.
#'
#' @return A modified `cell_plot` recipe.
#' @family cell-plot-modifiers
#'
#' @examples
#' se <- ReadPNA_Seurat(minimal_pna_pxl_file())
#' se <- LoadCellGraphs(se, cells = colnames(se)[4], verbose = FALSE) |>
#'   ComputeLayout(layout_method = "spectral")
#'
#' cell_graph <- CellGraphs(se)[[4]]
#'
#' layout_data <- FetchLayoutData(cell_graph, vars = "CD82", layout_method = "spectral_3d") |>
#'   # Downsample to speed up tests
#'   dplyr::slice_sample(n = 5000)
#'
#' cell_plot(layout_data) |>
#'   cell_theme(
#'     background_color = "navy",
#'     strip_background_color = "grey30",
#'     text_color = "white",
#'     text_size = 12
#'   )
#'
#' @export
cell_theme <- function(
  object,
  background_color,
  text_color,
  text_size,
  strip_background_color
) {
  defaults <- .cell_plot_default_theme
  if (missing(background_color)) {
    background_color <- defaults$background_color
  }
  if (missing(text_color)) {
    text_color <- defaults$text_color
  }
  if (missing(strip_background_color)) {
    strip_background_color <- defaults$strip_background_color
  }
  if (missing(text_size)) {
    text_size <- defaults$text_size
  }
  pixelatorR:::assert_valid_color(
    background_color,
    arg = "background_color"
  )
  pixelatorR:::assert_length(
    background_color,
    n = 1,
    arg_x = "background_color"
  )
  pixelatorR:::assert_valid_color(text_color, arg = "text_color")
  pixelatorR:::assert_length(text_color, n = 1, arg_x = "text_color")
  pixelatorR:::assert_valid_color(
    strip_background_color,
    arg = "strip_background_color"
  )
  pixelatorR:::assert_length(
    strip_background_color,
    n = 1,
    arg_x = "strip_background_color"
  )
  pixelatorR:::assert_single_value(text_size, type = "numeric", arg = "text_size")
  pixelatorR:::assert_within_limits(
    text_size,
    limits = c(.Machine$double.eps, Inf),
    arg = "text_size"
  )
  if (!is.finite(text_size)) {
    cli::cli_abort(c("i" = "{.arg text_size} must be finite."))
  }

  specification <- .cell_plot_default_theme
  specification$background_color <- background_color
  specification$strip_background_color <- strip_background_color
  specification$text_color <- text_color
  specification$text_size <- text_size
  return(.replace_cell_plot_field(object, "theme", specification))
}

#' Annotate a cell plot
#'
#' `r lifecycle::badge("experimental")`
#'
#' Sets title, subtitle, and color legend title text shared by cell plot
#' renderers. Annotation content is equivalent across backends, although
#' placement and typography may differ. When `legend_title` is omitted, the
#' mapped color column name is used.
#'
#' @param object A `cell_plot` recipe.
#' @param title,subtitle,legend_title Optional scalar character strings.
#'
#' @return A modified `cell_plot` recipe.
#' @family cell-plot-modifiers
#'
#' @examples
#' se <- ReadPNA_Seurat(minimal_pna_pxl_file())
#' se <- LoadCellGraphs(se, cells = colnames(se)[4], verbose = FALSE) |>
#'   ComputeLayout(layout_method = "spectral")
#'
#' cell_graph <- CellGraphs(se)[[4]]
#'
#' layout_data <- FetchLayoutData(cell_graph, vars = "CD82", layout_method = "spectral_3d") |>
#'   # Downsample to speed up tests
#'   dplyr::slice_sample(n = 5000)
#'
#' cell_plot(layout_data, color = CD82) |>
#'   cell_annotation(
#'     title = "CD82 distribution on a cell",
#'     subtitle = "Positive nodes colored in red",
#'     legend_title = "CD82 nodes"
#'   )
#'
#' @export
cell_annotation <- function(
  object,
  title = NULL,
  subtitle = NULL,
  legend_title = NULL
) {
  pixelatorR:::assert_single_value(
    title,
    type = "string",
    allow_null = TRUE,
    arg = "title"
  )
  pixelatorR:::assert_single_value(
    subtitle,
    type = "string",
    allow_null = TRUE,
    arg = "subtitle"
  )
  pixelatorR:::assert_single_value(
    legend_title,
    type = "string",
    allow_null = TRUE,
    arg = "legend_title"
  )
  if (is.null(title) && is.null(subtitle) && is.null(legend_title)) {
    cli::cli_abort(
      c(
        "i" = paste(
          "At least one of {.arg title}, {.arg subtitle}, or",
          "{.arg legend_title} must be supplied."
        )
      )
    )
  }

  specification <- .cell_plot_default_annotation
  specification["title"] <- list(title)
  specification["subtitle"] <- list(subtitle)
  specification["legend_title"] <- list(legend_title)
  return(.replace_cell_plot_field(object, "annotation", specification))
}

#' Add rotating coordinates to a cell plot
#'
#' `r lifecycle::badge("experimental")`
#'
#' Records rotation geometry for an animation renderer. This sets a rotation
#' sequence rather than a single static view. The number of frames belongs to
#' [cell_plot_animate()], not this modifier. Printed ggplot and Plotly views keep
#' the unrotated coordinates. The default axis is `"y"`, which spins the cell
#' left to right on screen.
#'
#' @param object A `cell_plot` recipe.
#' @param axis Axis around which coordinates are rotated. Use `"x"`, `"y"`, or
#' `"z"` for a principal axis, or three finite numeric values for an arbitrary
#' axis. Numeric axes are normalized before they are stored.
#' @param max_degree Maximum rotation angle from -360 to 360 degrees. Positive
#' angles follow the right-hand rule; use a negative angle to reverse direction.
#' @param boomerang Whether the animation returns through the frame sequence.
#' @param origin Rotation origin. `"origo"` uses `(0, 0, 0)`. `"centroid"`
#' calculates the centroid independently within each [cell_grid()] panel.
#'
#' @return A modified `cell_plot` recipe.
#'
#' @seealso [cell_plot_animate()], [cell_illuminate()]
#'
#' @examples
#' se <- ReadPNA_Seurat(minimal_pna_pxl_file())
#' se <- LoadCellGraphs(se, cells = colnames(se)[4], verbose = FALSE) |>
#'   ComputeLayout(layout_method = "spectral")
#'
#' cell_graph <- CellGraphs(se)[[4]]
#'
#' layout_data <- FetchLayoutData(cell_graph, vars = "CD82", layout_method = "spectral_3d") |>
#'   # Downsample to speed up tests
#'   dplyr::slice_sample(n = 5000)
#'
#' cell_plot(layout_data) |>
#'   cell_coord_rotate(axis = "y")
#'
#' @export
cell_coord_rotate <- function(
  object,
  axis = "y",
  max_degree = 360,
  boomerang = FALSE,
  origin = c("origo", "centroid")
) {
  axis <- if (is.character(axis)) {
    match.arg(axis, choices = c("x", "y", "z"))
  } else {
    .normalize_cell_direction(axis, argument = "axis")
  }
  pixelatorR:::assert_single_value(
    max_degree,
    type = "numeric",
    arg = "max_degree"
  )
  pixelatorR:::assert_within_limits(
    max_degree,
    limits = c(-360, 360),
    arg = "max_degree"
  )
  pixelatorR:::assert_single_value(
    boomerang,
    type = "bool",
    arg = "boomerang"
  )
  if (!is.finite(max_degree)) {
    cli::cli_abort(c("i" = "{.arg max_degree} must be finite."))
  }
  origin <- match.arg(origin)

  specification <- list(
    type = "rotate",
    axis = axis,
    max_degree = max_degree,
    boomerang = boomerang,
    origin = origin
  )
  return(.replace_cell_plot_field(object, "coord", specification))
}
