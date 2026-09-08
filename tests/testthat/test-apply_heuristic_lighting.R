test_that("heuristic_illumination works as expected", {
  set.seed(123)

  layout <- tibble::tibble(
    x = rnorm(100),
    y = rnorm(100),
    z = rnorm(100)
  )

  expect_no_error(
    illum <- heuristic_illumination(layout)
  )

  expect_type(illum, "double")
  expect_length(illum, nrow(layout))

  # deterministic regression check (allow tiny numeric drift)
  expect_equal(sum(illum), 56.83794, tolerance = 1e-4)

  expect_no_error(
    illum <- heuristic_illumination(layout,
      normalize_weights = FALSE
    )
  )
  expect_equal(sum(illum), 125.0435, tolerance = 1e-4)

  expect_no_error(
    illum <- heuristic_illumination(layout,
      directional_light_weight = 0.5,
      volume_shading_weight = 0.3,
      ambient_occlusion_weight = 0.2,
      normalize_weights = TRUE
    )
  )
  expect_equal(sum(illum), 48.59722, tolerance = 1e-4)

  # Default light_direction is +z
  expect_equal(
    heuristic_illumination(layout),
    heuristic_illumination(layout, light_direction = c(0, 0, 1))
  )
  expect_equal(
    heuristic_illumination(layout, light_direction = c(0, 0, 2)),
    heuristic_illumination(layout, light_direction = c(0, 0, 1))
  )

  # Opposite direction inverts the directional ranking when other weights are zero
  illum_pos <- heuristic_illumination(
    layout,
    volume_shading_weight = 0,
    ambient_occlusion_weight = 0,
    clamp_quantiles = c(0, 1),
    light_direction = c(0, 0, 1)
  )
  illum_neg <- heuristic_illumination(
    layout,
    volume_shading_weight = 0,
    ambient_occlusion_weight = 0,
    clamp_quantiles = c(0, 1),
    light_direction = c(0, 0, -1)
  )
  expect_equal(illum_neg, 1 - illum_pos)

  illum_diag <- heuristic_illumination(layout, light_direction = c(1, 1, 1))
  expect_true(all(is.finite(illum_diag)))
  expect_length(illum_diag, nrow(layout))

  # Expect errors with bad input

  expect_error(
    heuristic_illumination(
      layout,
      directional_light_weight = 0,
      volume_shading_weight = 0,
      ambient_occlusion_weight = 0,
      normalize_weights = TRUE
    ),
    "At least one weight must be > 0"
  )
  bad_missing <- tibble::tibble(x = rnorm(10), y = rnorm(10))
  expect_error(heuristic_illumination(bad_missing), "x|y|z|layout")

  bad_type <- tibble::tibble(
    x = rnorm(10),
    y = as.character(rnorm(10)),
    z = rnorm(10)
  )
  expect_error(heuristic_illumination(bad_type), "numeric|vector|type")

  infinite_x <- tibble::tibble(
    x = c(rnorm(9), Inf),
    y = rnorm(10),
    z = rnorm(10)
  )

  expect_error(
    heuristic_illumination(infinite_x),
    "finite values|finite"
  )

  expect_error(
    heuristic_illumination(layout, clamp_quantiles = c(0.8, 0.2)),
    "less than|clamp_quantiles"
  )

  expect_error(
    heuristic_illumination(layout, clamp_quantiles = c(-0.1, 0.9)),
    "0|1|limits|within"
  )

  expect_error(
    heuristic_illumination(layout, clamp_quantiles = c(0.1, 1.1)),
    "0|1|limits|within"
  )

  expect_error(
    heuristic_illumination(layout, directional_light_weight = -0.1),
    "0|Inf|non-negative|within"
  )
  expect_error(
    heuristic_illumination(layout, volume_shading_weight = -0.1),
    "0|Inf|non-negative|within"
  )
  expect_error(
    heuristic_illumination(layout, ambient_occlusion_weight = -0.1),
    "0|Inf|non-negative|within"
  )

  expect_error(
    heuristic_illumination(layout, ambient_occlusion_k = 0),
    "1|Inf|ambient_occlusion_k|within"
  )

  expect_error(
    heuristic_illumination(layout, light_direction = c(0, 0, 0)),
    "non-zero|zero"
  )
  expect_error(
    heuristic_illumination(layout, light_direction = c(1, 0)),
    "exactly|length|element"
  )
  expect_error(
    heuristic_illumination(layout, light_direction = c(1, 0, NA_real_)),
    "finite"
  )
  expect_error(
    heuristic_illumination(layout, light_direction = c(1, 0, Inf)),
    "finite"
  )
})
