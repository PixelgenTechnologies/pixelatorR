test_that("Cell plot animation geometry works as expected", {
  specification <- list(
    max_degree = 360,
    boomerang = FALSE
  )
  expect_equal(
    list(
      full_turn = .cell_animation_angles(specification, frames = 4),
      partial_turn = .cell_animation_angles(
        list(max_degree = -180, boomerang = FALSE),
        frames = 4
      ),
      even_boomerang = .cell_animation_angles(
        list(max_degree = 180, boomerang = TRUE),
        frames = 6
      ),
      odd_boomerang = .cell_animation_angles(
        list(max_degree = 180, boomerang = TRUE),
        frames = 5
      ),
      single_frame = .cell_animation_angles(specification, frames = 1)
    ),
    list(
      full_turn = c(0, 90, 180, 270),
      partial_turn = c(0, -60, -120, -180),
      even_boomerang = c(0, 60, 120, 180, 120, 60),
      odd_boomerang = c(0, 90, 180, 120, 60),
      single_frame = 0
    )
  )
  expect_error(.cell_animation_angles(specification, frames = 0))
  expect_error(.cell_animation_angles(specification, frames = 1.5))
  expect_equal(
    .cell_animation_plot_limits(
      list(x = c(-1, 1), y = c(0, 0))
    ),
    list(x = c(-1.1, 1.1), y = c(-0.05, 0.05))
  )

  layout <- tibble::tibble(
    x = c(1, 0),
    y = c(0, 1),
    z = c(0, 0)
  )
  expect_equal(
    .cell_rotate_layout(
      layout,
      axis = "z",
      angle = 90,
      origin = "origo",
      panel_rows = list(1:2)
    ),
    tibble::tibble(
      x = c(0, -1),
      y = c(1, 0),
      z = c(0, 0)
    ),
    tolerance = 1e-10
  )
  expect_equal(
    .cell_rotate_layout(
      layout[1, ],
      axis = c(1, 1, 1),
      angle = 120,
      origin = "origo",
      panel_rows = list(1L)
    ),
    tibble::tibble(x = 0, y = 1, z = 0),
    tolerance = 1e-10
  )

  panel_layout <- tibble::tibble(
    x = c(0, 2, 10, 14),
    y = 0,
    z = 0
  )
  expect_equal(
    .cell_rotate_layout(
      panel_layout,
      axis = "z",
      angle = 180,
      origin = "centroid",
      panel_rows = list(1:2, 3:4)
    ),
    tibble::tibble(
      x = c(2, 0, 14, 10),
      y = c(0, 0, 0, 0),
      z = c(0, 0, 0, 0)
    ),
    tolerance = 1e-10
  )

  plot_data <- tibble::tibble(
    id = c("left", "middle", "right"),
    x = c(-1, 0, 1),
    y = 0,
    z = 0,
    marker = c(10, 20, 30),
    confidence = c(0.1, 0.2, 0.3)
  )
  built <- cell_plot(
    plot_data,
    color = marker,
    alpha = confidence
  ) |>
    cell_node_scale_color(
      colors = c("black", "white"),
      limits = c(10, 30)
    ) |>
    cell_node_scale_alpha(
      alphas = c(0.2, 1),
      limits = c(0.1, 0.3)
    ) |>
    cell_coord_rotate(axis = "y") |>
    build_cell_plot()
  frame <- .cell_animation_frame(built, angle = 90)

  expect_equal(
    list(
      data = frame$data,
      color = frame$color[c("values", "resolved", "illuminated")],
      alpha = frame$alpha[c("values", "resolved")],
      depth_mapping = frame$mapping$depth,
      apparent_size = .cell_plot_projected_sizes(frame),
      limits = .cell_animation_limits(
        built,
        angles = c(0, 90, 180, 270)
      )
    ),
    list(
      data = tibble::tibble(
        id = c("right", "middle", "left"),
        x = c(0, 0, 0),
        y = c(0, 0, 0),
        z = c(-1, 0, 1),
        marker = c(30, 20, 10),
        confidence = c(0.3, 0.2, 0.1)
      ),
      color = list(
        values = c(30, 20, 10),
        resolved = c("#FFFFFF", "#777777", "#000000"),
        illuminated = NULL
      ),
      alpha = list(
        values = c(0.3, 0.2, 0.1),
        resolved = c(1, 0.6, 0.2)
      ),
      depth_mapping = "z",
      apparent_size = c(5 / 7, 1, 5 / 3),
      limits = list(x = c(-1, 1), y = c(0, 0))
    ),
    tolerance = 1e-10
  )
})

test_that("Animation frames preserve overlapping coordinate mappings", {
  plot_data <- tibble::tibble(
    shared = c(1, 0),
    y = c(0, 1)
  )
  built <- cell_plot(
    plot_data,
    x = shared,
    y = y,
    z = shared,
    depth = NULL
  ) |>
    cell_coord_rotate(axis = "y") |>
    build_cell_plot()

  frame <- .cell_animation_frame(built, angle = 90)
  coordinate_mappings <- unname(unlist(frame$mapping[c("x", "y", "z")]))

  expect_equal(
    list(
      mappings = coordinate_mappings,
      coordinates = frame$data[coordinate_mappings],
      depth = frame$mapping$depth
    ),
    list(
      mappings = c(
        ".cell_animation_x",
        ".cell_animation_y",
        ".cell_animation_z"
      ),
      coordinates = tibble::tibble(
        .cell_animation_x = c(1, 0),
        .cell_animation_y = c(0, 1),
        .cell_animation_z = c(-1, 0)
      ),
      depth = NULL
    ),
    tolerance = 1e-10
  )
})

test_that("Animation renderers and encoding work as expected", {
  plot_data <- tibble::tibble(
    x = c(-1, 0, 1),
    y = 0,
    z = 0,
    marker = c(0, 1, 2),
    panel = c("a", "a", "b")
  )
  recipe <- cell_plot(plot_data, color = marker) |>
    cell_grid(cols = panel) |>
    cell_node_scale_color(colors = c("black", "white"), limits = c(0, 2)) |>
    cell_annotation(title = "Cells", subtitle = "demo") |>
    cell_coord_rotate(axis = "y", max_degree = 90)

  expect_error(cell_plot_animate(cell_plot(plot_data), tempfile(fileext = ".gif")))
  expect_error(
    cell_plot_animate(
      recipe,
      file.path(tempdir(), "missing-dir", "cell.gif")
    )
  )

  built <- build_cell_plot(recipe)
  frame <- .cell_animation_frame(built, angle = 90)
  limits <- list(x = c(-2, 2), y = c(-2, 2))
  ggplot_frame <- .render_cell_plot_ggplot(frame, limits = limits)
  expect_s3_class(ggplot_frame, "ggplot")
  expect_equal(
    ggplot_frame$coordinates$limits,
    list(x = c(-2, 2), y = c(-2, 2))
  )

  png_file <- tempfile(fileext = ".png")
  grDevices::png(png_file, width = 240, height = 240)
  expect_no_error(.render_cell_plot_base(frame, limits = limits))
  grDevices::dev.off()
  expect_true(file.exists(png_file))
  unlink(png_file)
  expect_equal(
    .cell_base_layout_matrix(
      n_row = 2L,
      n_col = 2L,
      has_title = FALSE,
      has_subtitle = FALSE,
      has_row_strips = TRUE,
      has_col_strips = TRUE,
      has_legend = FALSE
    ),
    list(
      mat = matrix(
        c(1L, 4L, 7L, 2L, 5L, 8L, 3L, 6L, 9L),
        nrow = 3L
      ),
      widths = c(0.14, 1, 1),
      heights = c(0.14, 1, 1),
      respect = TRUE
    )
  )

  constant_color <- cell_plot(plot_data) |>
    cell_grid(cols = panel) |>
    cell_coord_rotate(axis = "y", max_degree = 90) |>
    build_cell_plot()
  constant_frame <- .cell_animation_frame(constant_color, angle = 90)
  constant_png <- tempfile(fileext = ".png")
  grDevices::png(constant_png, width = 240, height = 240)
  expect_no_error(.render_cell_plot_base(constant_frame, limits = limits))
  grDevices::dev.off()
  unlink(constant_png)

  png_dimensions <- function(path) {
    connection <- file(path, "rb")
    on.exit(close(connection), add = TRUE)
    header <- readBin(connection, what = "raw", n = 24)
    return(c(
      width = sum(as.integer(header[17:20]) * 256^(3:0)),
      height = sum(as.integer(header[21:24]) * 256^(3:0))
    ))
  }
  ggplot_png <- tempfile(fileext = ".png")
  .cell_animation_write_frame(
    object = built,
    angle = 90,
    path = ggplot_png,
    limits = limits,
    width = 160L,
    height = 160L,
    res = 72,
    frame_backend = "ggplot2",
    png_device = .cell_animation_png_device()
  )
  expect_equal(png_dimensions(ggplot_png), c(width = 160, height = 160))
  unlink(ggplot_png)

  aspect_file <- tempfile(fileext = ".png")
  grDevices::png(aspect_file, width = 300, height = 200)
  .cell_base_scatter(
    x = c(-1, 1),
    y = c(-1, 1),
    colors = c("black", "white"),
    cex = c(1, 1),
    limits = list(x = c(-2, 2), y = c(-1, 1)),
    background = "white"
  )
  user_coordinates <- graphics::par("usr")
  plot_size <- graphics::par("pin")
  grDevices::dev.off()
  expect_equal(
    diff(user_coordinates[1:2]) / plot_size[[1]],
    diff(user_coordinates[3:4]) / plot_size[[2]]
  )
  unlink(aspect_file)

  expect_equal(
    .cell_base_continuous_legend_scale(
      list(colors = c("black", "white"), limits = c(0, 2)),
      n = 3L
    ),
    list(
      fills = c("#000000", "#777777", "#FFFFFF"),
      labels = c(0, 2)
    )
  )

  skip_if_not_installed("gifski")
  gif_file <- tempfile(fileext = ".gif")
  expect_equal(
    cell_plot_animate(
      recipe,
      gif_file,
      frames = 2,
      width = 160,
      height = 160,
      res = 72,
      fps = 2
    ),
    gif_file
  )
  expect_true(file.exists(gif_file))
  expect_gt(file.size(gif_file), 0)
  expect_equal(
    cell_plot_animate(
      recipe,
      gif_file,
      frames = 2,
      width = 160,
      height = 160,
      res = 72,
      fps = 2,
      frame_backend = "ggplot2"
    ),
    gif_file
  )
  unlink(gif_file)

  skip_if_not_installed("av")
  mp4_file <- tempfile(fileext = ".mp4")
  expect_no_error(
    cell_plot_animate(
      recipe,
      mp4_file,
      frames = 2,
      width = 160,
      height = 160,
      res = 72,
      fps = 2,
      workers = 2L
    )
  )
  expect_true(file.exists(mp4_file))
  unlink(mp4_file)
})

test_that("Locked illumination follows rotated frame geometry", {
  skip_if_not_installed("FNN")

  plot_data <- tibble::tibble(
    x = c(-1, 0, 1),
    y = c(0, 1, 0),
    z = 0
  )
  locked <- cell_plot(plot_data) |>
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
  unlocked <- cell_plot(plot_data) |>
    cell_illuminate(
      clamp_quantiles = c(0, 1),
      directional_light_weight = 1,
      volume_shading_weight = 0,
      ambient_occlusion_weight = 0,
      light_direction = c(0, 0, 1)
    ) |>
    cell_coord_rotate(axis = "y") |>
    build_cell_plot()

  expect_equal(
    list(
      locked_baked = locked$color$illuminated,
      locked_frame = .cell_animation_frame(
        locked,
        angle = 90
      )$color$illuminated,
      unlocked_frame = .cell_animation_frame(
        unlocked,
        angle = 90
      )$color$illuminated
    ),
    list(
      locked_baked = c("#333333", "#333333", "#333333"),
      locked_frame = c("#0F0F0F", "#212121", "#333333"),
      unlocked_frame = c("#333333", "#333333", "#333333")
    )
  )
})
