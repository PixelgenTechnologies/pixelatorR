test_that("Cell illumination components match expected results", {
  skip_if_not_installed("FNN")

  layout <- tibble::tibble(
    x = c(0, 1, 0, -1),
    y = c(0, 0, 1, 0),
    z = c(-1, 0, 1, 0.5)
  )
  coordinates <- .cell_illumination_coordinates(layout)
  components <- list(
    directional = .cell_directional_illumination(coordinates),
    volume = .cell_volume_illumination(coordinates),
    ambient = .cell_ambient_occlusion_illumination(coordinates, k = 2L)
  )
  combined <- .cell_combine_illumination(
    directional = components$directional,
    volume = components$volume,
    ambient = components$ambient,
    weights = c(0.7, 0.5, 1),
    clamp_quantiles = c(0, 1),
    normalize_weights = TRUE
  )
  local_illumination <- .cell_heuristic_illumination(
    layout,
    clamp_quantiles = c(0, 1),
    ambient_occlusion_k = 2L
  )

  expect_equal(
    list(
      components = components,
      combined = combined,
      local = local_illumination,
      legacy = pixelatorR::heuristic_illumination(
        layout,
        clamp_quantiles = c(0, 1),
        ambient_occlusion_k = 2L
      )
    ),
    list(
      components = list(
        directional = c(0, 0.5, 1, 0.75),
        volume = c(0, 0, 1, 0.28495925646099),
        ambient = c(
          0.571911575961696,
          1,
          0.428088424038312,
          0
        )
      ),
      combined = c(
        0.259959807255316,
        0.613636363636364,
        0.740040192744687,
        0.303399831013861
      ),
      local = c(
        0.259959807255316,
        0.613636363636364,
        0.740040192744687,
        0.303399831013861
      ),
      legacy = c(
        0.259959807255316,
        0.613636363636364,
        0.740040192744687,
        0.303399831013861
      )
    ),
    tolerance = 1e-12
  )
})

test_that("Cell illumination applies cached color channels", {
  expect_equal(
    list(
      hsv = .cell_apply_hsv_illumination(
        c("#6699CC", "#FF0000"),
        illumination = c(0.2, 0.8),
        ambient_intensity = 0.3,
        saturation_boost = 0.6
      ),
      hsv_cached = .cell_apply_hsv_channels(
        .cell_hsv_channels(c("#6699CC", "#FF0000")),
        illumination = c(0.2, 0.8),
        ambient_intensity = 0.3,
        saturation_boost = 0.6
      ),
      palette = .cell_apply_palette_illumination(
        c("#6699CC", "#FF0000"),
        shadow_colors = c("#000000", "#FFFFFF"),
        illumination = c(0.2, 0.8),
        ambient_intensity = 0.3
      )
    ),
    list(
      hsv = c("#1E3C5A", "#DB0000"),
      hsv_cached = c("#1E3C5A", "#DB0000"),
      palette = c("#2C4359", "#FF2323")
    )
  )
})

test_that("Illumination masks preserve excluded node colors", {
  skip_if_not_installed("FNN")

  plot_data <- tibble::tibble(
    x = c(-1, 0, 1),
    y = c(0, 1, 0),
    z = c(-1, 0, 1),
    detected = c(FALSE, TRUE, TRUE)
  )
  built <- cell_plot(plot_data, illumination_mask = detected) |>
    cell_illuminate(
      clamp_quantiles = c(0, 1),
      directional_light_weight = 1,
      volume_shading_weight = 0,
      ambient_occlusion_weight = 0,
      light_direction = c(0, 0, 1),
      lock_light = TRUE
    ) |>
    cell_coord_rotate(axis = "y") |>
    build_cell_plot()
  prepared <- .cell_prepare_animation_illumination(built)
  frame <- .cell_animation_frame(built, angle = -90)

  expect_equal(
    list(
      mapping = built$mapping$illumination_mask,
      resolved = built$color$resolved,
      illuminated = built$color$illuminated,
      cached_positions = attr(
        prepared,
        "cell_illumination_cache"
      )$panels[[1]]$apply_positions,
      frame_detected = frame$data$detected,
      frame_colors = frame$color$illuminated
    ),
    list(
      mapping = "detected",
      resolved = "#E5E5E5",
      illuminated = c("#E5E5E5", "#959595", "#E5E5E5"),
      cached_positions = c(2L, 3L),
      frame_detected = c(FALSE, TRUE, TRUE),
      frame_colors = c("#E5E5E5", "#959595", "#E5E5E5")
    )
  )
})

test_that("Illumination masks survive dropping missing size rows", {
  plot_data <- tibble::tibble(
    x = c(-1, 0, 1, 2),
    y = c(0, 1, 0, 0),
    z = c(-1, 10, 0, 1),
    size = c(1, NA, 1, 1),
    detected = c(FALSE, TRUE, TRUE, TRUE)
  )

  expect_warning(
    built <- cell_plot(
      plot_data,
      size = size,
      illumination_mask = detected
    ) |>
      cell_illuminate(
        clamp_quantiles = c(0, 1),
        directional_light_weight = 1,
        volume_shading_weight = 0,
        ambient_occlusion_weight = 0,
        light_direction = c(0, 0, 1)
      ) |>
      build_cell_plot(),
    "Removed 1 row containing missing values in size and/or alpha."
  )

  expect_equal(
    list(
      z = built$data$z,
      detected = built$data$detected,
      resolved = built$color$resolved,
      illuminated = built$color$illuminated
    ),
    list(
      z = c(-1, 0, 1),
      detected = c(FALSE, TRUE, TRUE),
      resolved = "#E5E5E5",
      illuminated = c("#E5E5E5", "#959595", "#E5E5E5")
    )
  )
})

test_that("Locked illumination reuses rotation-invariant terms", {
  skip_if_not_installed("FNN")

  plot_data <- tibble::tibble(
    x = c(-1, 0, 1, 0),
    y = c(0, 1, 0, -1),
    z = c(-0.5, 0, 0.5, 1)
  )
  make_built <- function(origin) {
    return(
      cell_plot(plot_data) |>
        cell_illuminate(
          clamp_quantiles = c(0, 1),
          ambient_occlusion_k = 2L,
          lock_light = TRUE
        ) |>
        cell_coord_rotate(axis = "y", origin = origin) |>
        build_cell_plot() |>
        .cell_prepare_animation_illumination()
    )
  }
  origo <- make_built("origo")
  centroid <- make_built("centroid")
  origo_cache <- attr(origo, "cell_illumination_cache")$panels[[1]]
  centroid_cache <- attr(centroid, "cell_illumination_cache")$panels[[1]]

  expected_frame_colors <- function(object, angle) {
    mapping <- object$mapping
    panel_rows <- .cell_plot_panel_rows(object$data, object$grid)
    layout <- tibble::tibble(
      x = object$data[[mapping$x]],
      y = object$data[[mapping$y]],
      z = object$data[[mapping$z]]
    )
    rotated <- .cell_rotate_layout(
      layout = layout,
      axis = object$coord$axis,
      angle = angle,
      origin = object$coord$origin,
      panel_rows = panel_rows
    )
    colors <- .cell_illuminated_colors(
      layout = rotated,
      colors = object$color$resolved,
      mapped_values = object$color$values,
      panel_rows = panel_rows,
      specification = object$illuminate
    )
    return(colors[order(rotated$z, na.last = TRUE)])
  }

  expect_equal(
    list(
      origo_volume_cached = !is.null(origo_cache$volume),
      centroid_volume_cached = !is.null(centroid_cache$volume),
      ambient_rows = length(origo_cache$ambient),
      color_channel_dimensions = dim(origo_cache$color_channels),
      origo_frame = .cell_animation_frame(
        origo,
        angle = 90
      )$color$illuminated,
      origo_expected = expected_frame_colors(origo, angle = 90),
      centroid_frame = .cell_animation_frame(
        centroid,
        angle = 90
      )$color$illuminated,
      centroid_expected = expected_frame_colors(centroid, angle = 90)
    ),
    list(
      origo_volume_cached = TRUE,
      centroid_volume_cached = FALSE,
      ambient_rows = 4L,
      color_channel_dimensions = c(3L, 4L),
      origo_frame = expected_frame_colors(origo, angle = 90),
      origo_expected = expected_frame_colors(origo, angle = 90),
      centroid_frame = expected_frame_colors(centroid, angle = 90),
      centroid_expected = expected_frame_colors(centroid, angle = 90)
    )
  )
})
