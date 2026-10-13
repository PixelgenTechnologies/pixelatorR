#' Rescale an illumination term safely
#'
#' Constant or non-finite terms are replaced by the midpoint of the requested
#' output range.
#'
#' @param values Numeric illumination values.
#' @param to Numeric output range.
#'
#' @return A numeric vector with the same length as `values`.
#'
#' @noRd
.cell_safe_rescale <- function(values, to = c(0, 1)) {
  value_range <- range(values, na.rm = TRUE)
  if (
    !all(is.finite(value_range)) ||
      diff(value_range) == 0
  ) {
    return(rep(mean(to), length(values)))
  }
  return(scales::rescale(values, to = to, from = value_range))
}

#' Extract cell illumination coordinates
#'
#' @param layout A data frame with numeric `x`, `y`, and `z` columns.
#'
#' @return A finite numeric matrix with three columns.
#'
#' @noRd
.cell_illumination_coordinates <- function(layout) {
  assert_class(layout, c("data.frame", "tbl_df"))
  assert_x_in_y(x = c("x", "y", "z"), y = names(layout))
  for (coordinate in c("x", "y", "z")) {
    assert_vector(layout[[coordinate]], type = "numeric")
  }
  coordinates <- as.matrix(layout[, c("x", "y", "z")])
  if (any(!is.finite(coordinates))) {
    cli::cli_abort("Columns {.field x}, {.field y}, and {.field z} must contain only finite values.")
  }
  return(coordinates)
}

#' Calculate directional cell illumination
#'
#' @param coordinates A numeric matrix with x, y, and z columns.
#' @param light_direction A three-dimensional light direction.
#'
#' @return Directional illumination scaled between zero and one.
#'
#' @noRd
.cell_directional_illumination <- function(
  coordinates,
  light_direction = c(0, 0, 1)
) {
  light_direction <- .normalize_cell_direction(
    light_direction,
    argument = "light_direction"
  )
  directional <- as.numeric(coordinates %*% light_direction)
  return(.cell_safe_rescale(directional))
}

#' Calculate radial cell illumination
#'
#' @param coordinates A numeric matrix with x, y, and z columns.
#'
#' @return Radial illumination scaled between zero and one.
#'
#' @noRd
.cell_volume_illumination <- function(coordinates) {
  radius <- sqrt(rowSums(coordinates^2))
  return(.cell_safe_rescale(radius))
}

#' Calculate ambient-occlusion cell illumination
#'
#' @param coordinates A numeric matrix with x, y, and z columns.
#' @param k Positive number of nearest neighbors.
#'
#' @return Ambient occlusion scaled from one to zero.
#'
#' @noRd
.cell_ambient_occlusion_illumination <- function(coordinates, k) {
  expect_FNN()
  neighbor_distances <- FNN::get.knn(coordinates, k = k)$nn.dist
  ambient_occlusion <- rowMeans(sqrt(neighbor_distances))
  return(.cell_safe_rescale(ambient_occlusion, to = c(1, 0)))
}

#' Combine independent cell illumination terms
#'
#' @param directional,volume,ambient Numeric illumination terms of equal
#' length.
#' @param weights Three non-negative weights in directional, volume, and
#' ambient order.
#' @param clamp_quantiles Two ordered quantiles between zero and one.
#' @param normalize_weights Whether to scale weights to sum to one.
#'
#' @return A weighted and quantile-clamped illumination vector.
#'
#' @noRd
.cell_combine_illumination <- function(
  directional,
  volume,
  ambient,
  weights,
  clamp_quantiles = c(0.01, 0.95),
  normalize_weights = TRUE
) {
  term_lengths <- lengths(list(directional, volume, ambient))
  if (length(unique(term_lengths)) != 1L) {
    cli::cli_abort("Illumination terms must have equal lengths.")
  }
  assert_vector(weights, type = "numeric", n = 3)
  assert_within_limits(weights, limits = c(0, Inf))
  assert_vector(clamp_quantiles, type = "numeric", n = 2)
  assert_within_limits(clamp_quantiles, limits = c(0, 1))
  if (clamp_quantiles[[1]] >= clamp_quantiles[[2]]) {
    cli::cli_abort("{.arg clamp_quantiles[1]} must be less than {.arg clamp_quantiles[2]}.")
  }
  assert_single_value(
    normalize_weights,
    type = "bool",
    arg = "normalize_weights"
  )

  if (normalize_weights) {
    weight_sum <- sum(weights)
    if (weight_sum == 0) {
      cli::cli_abort("At least one illumination weight must be positive.")
    }
    weights <- weights / weight_sum
  }

  illumination <- weights[[1]] * directional +
    weights[[2]] * volume +
    weights[[3]] * ambient
  clamp <- stats::quantile(
    illumination,
    probs = clamp_quantiles,
    na.rm = TRUE
  )
  illumination <- pmin(
    pmax(illumination, clamp[[1]]),
    clamp[[2]]
  )
  return(illumination)
}

#' Calculate heuristic cell illumination
#'
#' Calculates directional, radial-volume, and ambient-occlusion terms
#' independently, then combines them using the requested weights.
#'
#' @param layout A data frame with numeric `x`, `y`, and `z` columns.
#' @param clamp_quantiles Two ordered clamp quantiles.
#' @param directional_light_weight,volume_shading_weight,ambient_occlusion_weight
#' Non-negative illumination weights.
#' @param ambient_occlusion_k Positive number of nearest neighbors.
#' @param normalize_weights Whether to scale weights to sum to one.
#' @param light_direction A three-dimensional light direction. The neutral
#' default mirrors [pixelatorR::heuristic_illumination()] so the two
#' implementations can be compared directly. The default lamp angle used by
#' cell plots belongs to [cell_illuminate()].
#'
#' @return A numeric illumination vector.
#'
#' @noRd
.cell_heuristic_illumination <- function(
  layout,
  clamp_quantiles = c(0.01, 0.95),
  directional_light_weight = 0.7,
  volume_shading_weight = 0.5,
  ambient_occlusion_weight = 1,
  ambient_occlusion_k = 20,
  normalize_weights = TRUE,
  light_direction = c(0, 0, 1)
) {
  coordinates <- .cell_illumination_coordinates(layout)
  rows <- nrow(coordinates)
  assert_single_value(
    ambient_occlusion_k,
    type = "integer",
    arg = "ambient_occlusion_k"
  )
  assert_within_limits(
    ambient_occlusion_k,
    limits = c(1, rows - 1L),
    arg = "ambient_occlusion_k"
  )
  weights <- c(
    directional_light_weight,
    volume_shading_weight,
    ambient_occlusion_weight
  )
  assert_vector(weights, type = "numeric", n = 3)
  assert_within_limits(weights, limits = c(0, Inf))
  zero_term <- numeric(rows)
  directional <- if (directional_light_weight > 0) {
    .cell_directional_illumination(coordinates, light_direction)
  } else {
    zero_term
  }
  volume <- if (volume_shading_weight > 0) {
    .cell_volume_illumination(coordinates)
  } else {
    zero_term
  }
  ambient <- if (ambient_occlusion_weight > 0) {
    .cell_ambient_occlusion_illumination(
      coordinates,
      k = ambient_occlusion_k
    )
  } else {
    zero_term
  }
  return(.cell_combine_illumination(
    directional = directional,
    volume = volume,
    ambient = ambient,
    weights = weights,
    clamp_quantiles = clamp_quantiles,
    normalize_weights = normalize_weights
  ))
}

#' Normalize a combined illumination mask
#'
#' @param illumination Numeric illumination values.
#'
#' @return Values scaled between zero and one. Constant values become one.
#'
#' @noRd
.cell_normalize_illumination <- function(illumination) {
  illumination_range <- range(illumination)
  if (diff(illumination_range) == 0) {
    return(rep(1, length(illumination)))
  }
  return(scales::rescale(
    illumination,
    to = c(0, 1),
    from = illumination_range
  ))
}

#' Parse cell colors into HSV channels
#'
#' @param hex_colors Valid hexadecimal colors.
#'
#' @return A three-row matrix containing hue, saturation, and value.
#'
#' @noRd
.cell_hsv_channels <- function(hex_colors) {
  return(grDevices::rgb2hsv(grDevices::col2rgb(hex_colors)))
}

#' Apply illumination to cached HSV channels
#'
#' @param hsv_channels A three-row hue, saturation, and value matrix.
#' @param illumination Numeric mask between zero and one.
#' @param ambient_intensity Minimum brightness.
#' @param saturation_boost Saturation increase in shadow.
#'
#' @return Illuminated hexadecimal colors.
#'
#' @noRd
.cell_apply_hsv_channels <- function(
  hsv_channels,
  illumination,
  ambient_intensity = 0.1,
  saturation_boost = 0.7
) {
  light_factor <- ambient_intensity +
    illumination * (1 - ambient_intensity)
  value <- hsv_channels[3, ] * light_factor
  saturation <- hsv_channels[2, ] *
    (1 + (1 - light_factor) * saturation_boost)
  return(grDevices::hsv(
    h = hsv_channels[1, ],
    s = pmin(1, saturation),
    v = value
  ))
}

#' Apply HSV illumination to cell colors
#'
#' Colors are already validated and resolved by the cell color scale, so this
#' helper deliberately avoids validating every color again.
#'
#' @param hex_colors Valid hexadecimal colors.
#' @param illumination Numeric mask between zero and one.
#' @param ambient_intensity Minimum brightness.
#' @param saturation_boost Saturation increase in shadow.
#'
#' @return Illuminated hexadecimal colors.
#'
#' @noRd
.cell_apply_hsv_illumination <- function(
  hex_colors,
  illumination,
  ambient_intensity = 0.1,
  saturation_boost = 0.7
) {
  hsv_channels <- .cell_hsv_channels(hex_colors)
  return(.cell_apply_hsv_channels(
    hsv_channels = hsv_channels,
    illumination = illumination,
    ambient_intensity = ambient_intensity,
    saturation_boost = saturation_boost
  ))
}

#' Parse cell colors into RGB channels
#'
#' @param hex_colors Valid hexadecimal colors.
#'
#' @return A three-row red, green, and blue matrix.
#'
#' @noRd
.cell_rgb_channels <- function(hex_colors) {
  return(grDevices::col2rgb(hex_colors))
}

#' Apply palette illumination to cached RGB channels
#'
#' @param rgb_channels A three-row red, green, and blue matrix.
#' @param shadow_colors Hexadecimal shadow colors.
#' @param illumination Numeric mask between zero and one.
#' @param ambient_intensity Minimum brightness.
#'
#' @return Illuminated hexadecimal colors.
#'
#' @noRd
.cell_apply_palette_channels <- function(
  rgb_channels,
  shadow_colors,
  illumination,
  ambient_intensity = 0.1
) {
  if (length(shadow_colors) == 1L) {
    shadow_channels <- matrix(
      grDevices::col2rgb(shadow_colors),
      nrow = 3,
      ncol = ncol(rgb_channels)
    )
  } else {
    shadow_channels <- grDevices::col2rgb(shadow_colors)
  }
  light_factor <- ambient_intensity +
    illumination * (1 - ambient_intensity)
  light_matrix <- matrix(
    light_factor,
    nrow = 3,
    ncol = length(light_factor),
    byrow = TRUE
  )
  blended <- rgb_channels * light_matrix +
    shadow_channels * (1 - light_matrix)
  return(grDevices::rgb(
    red = blended[1, ],
    green = blended[2, ],
    blue = blended[3, ],
    maxColorValue = 255
  ))
}

#' Apply palette illumination to cell colors
#'
#' Colors are already validated and resolved by the cell color scale, so this
#' helper deliberately avoids validating every color again.
#'
#' @param base_colors Valid hexadecimal base colors.
#' @param shadow_colors Valid hexadecimal shadow colors.
#' @param illumination Numeric mask between zero and one.
#' @param ambient_intensity Minimum brightness.
#'
#' @return Illuminated hexadecimal colors.
#'
#' @noRd
.cell_apply_palette_illumination <- function(
  base_colors,
  shadow_colors,
  illumination,
  ambient_intensity = 0.1
) {
  rgb_channels <- .cell_rgb_channels(base_colors)
  return(.cell_apply_palette_channels(
    rgb_channels = rgb_channels,
    shadow_colors = shadow_colors,
    illumination = illumination,
    ambient_intensity = ambient_intensity
  ))
}

#' Prepare reusable illumination terms for rotating frames
#'
#' Ambient occlusion is invariant under rigid rotation. Radial volume shading
#' is also invariant when rotation is around the coordinate origin. Base color
#' channels are parsed once for reuse by every frame.
#'
#' @param object A built cell plot.
#'
#' @return `object` with an internal illumination-cache attribute.
#'
#' @noRd
.cell_prepare_animation_illumination <- function(object) {
  if (
    is.null(object$illuminate) ||
      !isTRUE(object$illuminate$lock_light) ||
      !is.null(attr(object, "cell_illumination_cache", exact = TRUE))
  ) {
    return(object)
  }

  mapping <- object$mapping
  layout <- tibble::tibble(
    x = object$data[[mapping$x]],
    y = object$data[[mapping$y]],
    z = object$data[[mapping$z]]
  )
  colors <- rep_len(object$color$resolved, nrow(layout))
  panel_rows <- .cell_plot_panel_rows(object$data, object$grid)
  specification <- object$illuminate

  panels <- lapply(panel_rows, function(rows) {
    if (length(rows) < 2L) {
      return(list(rows = rows, skip = TRUE))
    }

    coordinates <- .cell_illumination_coordinates(
      layout[rows, , drop = FALSE]
    )
    ambient <- if (specification$ambient_occlusion_weight > 0) {
      .cell_ambient_occlusion_illumination(
        coordinates,
        k = min(specification$ambient_occlusion_k, length(rows) - 1L)
      )
    } else {
      numeric(length(rows))
    }
    volume <- if (
      specification$volume_shading_weight > 0 &&
        identical(object$coord$origin, "origo")
    ) {
      .cell_volume_illumination(coordinates)
    } else {
      NULL
    }
    apply_positions <- if (is.null(object$color$values)) {
      seq_along(rows)
    } else {
      which(!is.na(object$color$values[rows]))
    }
    if (!is.null(mapping$illumination_mask)) {
      apply_positions <- intersect(
        apply_positions,
        which(object$data[[mapping$illumination_mask]][rows])
      )
    }
    if (length(apply_positions) == 0L) {
      return(list(rows = rows, skip = TRUE))
    }
    base_colors <- colors[rows[apply_positions]]
    color_channels <- if (is.null(specification$shadow_colors)) {
      .cell_hsv_channels(base_colors)
    } else {
      .cell_rgb_channels(base_colors)
    }

    return(list(
      rows = rows,
      skip = FALSE,
      ambient = ambient,
      volume = volume,
      apply_positions = apply_positions,
      color_channels = color_channels
    ))
  })

  attr(object, "cell_illumination_cache") <- list(panels = panels)
  return(object)
}

#' Apply cached illumination to a rotated frame
#'
#' @param layout Rotated x, y, and z coordinates.
#' @param colors Unilluminated hexadecimal colors.
#' @param specification An illumination specification.
#' @param cache Reusable per-panel terms and color channels.
#'
#' @return One illuminated hexadecimal color per layout row.
#'
#' @noRd
.cell_apply_cached_illumination <- function(
  layout,
  colors,
  specification,
  cache
) {
  illuminated_colors <- rep_len(colors, nrow(layout))
  weights <- c(
    specification$directional_light_weight,
    specification$volume_shading_weight,
    specification$ambient_occlusion_weight
  )

  for (panel in cache$panels) {
    if (isTRUE(panel$skip)) {
      next
    }

    coordinates <- .cell_illumination_coordinates(
      layout[panel$rows, , drop = FALSE]
    )
    zero_term <- numeric(length(panel$rows))
    directional <- if (specification$directional_light_weight > 0) {
      .cell_directional_illumination(
        coordinates,
        specification$light_direction
      )
    } else {
      zero_term
    }
    volume <- if (specification$volume_shading_weight == 0) {
      zero_term
    } else if (!is.null(panel$volume)) {
      panel$volume
    } else {
      .cell_volume_illumination(coordinates)
    }
    illumination <- .cell_combine_illumination(
      directional = directional,
      volume = volume,
      ambient = panel$ambient,
      weights = weights,
      clamp_quantiles = specification$clamp_quantiles,
      normalize_weights = TRUE
    )
    illumination <- .cell_normalize_illumination(illumination)
    apply_mask <- illumination[panel$apply_positions]
    color_rows <- panel$rows[panel$apply_positions]

    if (is.null(specification$shadow_colors)) {
      illuminated_colors[color_rows] <- .cell_apply_hsv_channels(
        hsv_channels = panel$color_channels,
        illumination = apply_mask,
        ambient_intensity = specification$ambient_intensity,
        saturation_boost = specification$saturation_boost
      )
    } else {
      shadow_colors <- scales::col_numeric(
        palette = specification$shadow_colors,
        domain = c(0, 1),
        na.color = "transparent"
      )(apply_mask)
      illuminated_colors[color_rows] <- .cell_apply_palette_channels(
        rgb_channels = panel$color_channels,
        shadow_colors = shadow_colors,
        illumination = apply_mask,
        ambient_intensity = specification$ambient_intensity
      )
    }
  }

  return(illuminated_colors)
}
