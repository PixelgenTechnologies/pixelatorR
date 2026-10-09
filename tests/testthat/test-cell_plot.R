test_that("Cell plot recipes and modifiers work as expected", {
  plot_data <- tibble::tibble(
    x = c(0, 1),
    y = c(1, 0),
    z = c(-1, 1),
    marker = c(0, 2),
    abundance = c(1, 5),
    confidence = c(0.2, 0.8),
    cell = c("a", "b")
  )

  base_plot <- cell_plot(
    plot_data,
    color = marker,
    size = "abundance",
    alpha = confidence
  )

  expect_equal(
    base_plot,
    structure(
      list(
        data = plot_data,
        mapping = list(
          x = "x",
          y = "y",
          z = "z",
          depth = "z",
          color = "marker",
          size = "abundance",
          alpha = "confidence",
          illumination_mask = NULL,
          arrange = "z"
        ),
        constant = list(color = NULL, size = NULL, alpha = NULL),
        grid = NULL,
        color = NULL,
        size = NULL,
        depth = NULL,
        alpha = NULL,
        theme = NULL,
        annotation = NULL,
        illuminate = NULL,
        coord = NULL
      ),
      class = "cell_plot"
    )
  )
  expect_invisible(.validate_cell_plot(base_plot))

  invalid_plot <- base_plot
  invalid_plot$grid <- "cell"
  expect_error(.validate_cell_plot(invalid_plot))

  modified_plot <- base_plot |>
    cell_grid(cols = cell) |>
    cell_node_scale_color(
      colors = c("blue", "red"),
      limits = c(0, 2)
    ) |>
    cell_node_scale_size(sizes = c(2, 6), limits = c(1, 5)) |>
    cell_node_scale_alpha(alphas = c(0.2, 1), limits = c(0, 1)) |>
    cell_illuminate(ambient_occlusion_k = 10) |>
    cell_theme(
      background_color = "navy",
      strip_background_color = "grey30",
      text_color = "white",
      text_size = 12
    ) |>
    cell_annotation(
      title = "CD3 polarized cells",
      subtitle = "Selected cells",
      legend_title = "CD3"
    ) |>
    cell_coord_rotate(
      axis = "z",
      max_degree = 180
    )

  expect_equal(
    unclass(modified_plot)[-c(1, 2, 3)],
    list(
      grid = list(rows = NULL, cols = "cell"),
      color = list(
        colors = c("blue", "red"),
        limits = c(0, 2),
        na_color = "grey50",
        type = "auto"
      ),
      size = list(sizes = c(2, 6), limits = c(1, 5)),
      depth = NULL,
      alpha = list(alphas = c(0.2, 1), limits = c(0, 1)),
      theme = list(
        background_color = "navy",
        strip_background_color = "grey30",
        text_color = "white",
        text_size = 12
      ),
      annotation = list(
        title = "CD3 polarized cells",
        subtitle = "Selected cells",
        legend_title = "CD3"
      ),
      illuminate = list(
        clamp_quantiles = c(0.01, 0.95),
        directional_light_weight = 0.7,
        volume_shading_weight = 0.5,
        ambient_occlusion_weight = 1,
        ambient_occlusion_k = 10,
        ambient_intensity = 0.3,
        saturation_boost = 0.6,
        shadow_colors = NULL,
        light_direction = c(
          -0.639602149066831,
          0.426401432711221,
          0.639602149066831
        ),
        lock_light = FALSE
      ),
      coord = list(
        type = "rotate",
        axis = "z",
        max_degree = 180,
        origin = "origo"
      )
    )
  )
  expect_null(base_plot$grid)
  expect_null(base_plot$theme)

  last_theme_wins <- modified_plot |>
    cell_theme(
      background_color = "white",
      strip_background_color = "grey80",
      text_color = "black",
      text_size = 10
    )
  expect_equal(
    last_theme_wins$theme,
    list(
      background_color = "white",
      strip_background_color = "grey80",
      text_color = "black",
      text_size = 10
    )
  )

  positional_theme <- cell_theme(base_plot, "navy", "white", 12)
  expect_equal(
    positional_theme$theme,
    list(
      background_color = "navy",
      strip_background_color = "#D9D9D9",
      text_color = "white",
      text_size = 12
    )
  )

  expect_error(cell_plot(dplyr::select(plot_data, -z), color = "marker"))
  expect_error(cell_plot(plot_data, z = NULL))

  plot_summary <- summary(base_plot)
  expected_summary <- c(
    "<cell_plot>",
    "Rows: 2",
    "Mappings: x, y, z, depth, color, size, alpha, arrange",
    "Constants: none",
    "Modifiers: none"
  )
  expect_s3_class(plot_summary, "summary.cell_plot")
  expect_equal(as.character(plot_summary), expected_summary)
  expect_equal(capture.output(print(plot_summary)), expected_summary)
  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  expect_s3_class(print(base_plot), "ggplot")
  grDevices::dev.off()
  unlink(plot_file)

  expect_error(cell_plot(as.data.frame(plot_data)))
  expect_error(cell_plot(plot_data[0, ]))
  expect_error(cell_plot(plot_data, color = log(marker)))
  expect_error(
    cell_plot(dplyr::mutate(plot_data, marker = NA_real_), color = marker)
  )
  expect_error(
    cell_plot(dplyr::mutate(plot_data, abundance = NA_character_), size = abundance)
  )
  expect_error(
    cell_plot(dplyr::mutate(plot_data, confidence = Inf), alpha = confidence),
    "finite"
  )
  expect_error(
    cell_plot(dplyr::mutate(plot_data, marker = c(1, Inf)), color = marker),
    "finite"
  )
  expect_error(cell_plot(plot_data, illumination_mask = marker))
  expect_error(
    cell_plot(
      dplyr::mutate(plot_data, detected = c(TRUE, NA)),
      illumination_mask = detected
    )
  )
  expect_error(cell_grid(base_plot))
  expect_error(cell_grid(base_plot, cols = marker))
  grid_only <- cell_grid(base_plot, rows = cell, cols = cell)
  expect_equal(nrow(grid_only$data), nrow(base_plot$data))
  expect_equal(grid_only$data, base_plot$data)
  expect_equal(
    pixelatorR:::.cell_plot_panel_rows(plot_data, grid_only$grid),
    structure(
      list(1L, 2L),
      ptype = integer(0),
      class = c("vctrs_list_of", "vctrs_vctr", "list")
    )
  )
  expect_error(
    cell_node_scale_color(base_plot, colors = "not-a-color")
  )
  expect_error(
    cell_node_scale_color(
      base_plot,
      colors = c("blue", "red"),
      limits = c(2, 0)
    )
  )
  expect_error(
    cell_plot(plot_data, color = cell) |>
      cell_node_scale_color(colors = c(a = "blue"))
  )
  expect_error(cell_node_scale_size(base_plot, sizes = c(6, 2)))
  expect_error(cell_node_scale_size(base_plot, sizes = c(2, 6), limits = c(5, 1)))
  expect_error(cell_node_scale_alpha(base_plot, alphas = c(-0.1, 1)))
  expect_error(cell_node_scale_alpha(base_plot, alphas = c(0.2, 1), limits = c(1, 0)))
  expect_error(cell_illuminate(base_plot, ambient_intensity = 1.5))
  expect_error(cell_illuminate(base_plot, saturation_boost = Inf))
  expect_error(cell_illuminate(base_plot, shadow_colors = "not-a-color"))
  expect_error(cell_illuminate(base_plot, light_direction = c(0, 0, 0)))
  expect_error(cell_illuminate(base_plot, light_direction = c(0, 1)))
  expect_error(cell_illuminate(base_plot, lock_light = "yes"))
  expect_error(
    cell_illuminate(
      base_plot,
      directional_light_weight = 0,
      volume_shading_weight = 0,
      ambient_occlusion_weight = 0
    ),
    "positive"
  )
  expect_error(
    cell_illuminate(base_plot, clamp_quantiles = c(0.9, 0.1)),
    "increasing"
  )
  expect_error(cell_illuminate(base_plot, ambient_occlusion_k = 0))
  expect_error(cell_theme(base_plot, text_size = 0))
  expect_error(cell_theme(base_plot, strip_background_color = "not-a-color"))
  expect_error(cell_annotation(base_plot))
  expect_equal(
    cell_annotation(base_plot, legend_title = "CD3")$annotation,
    list(
      title = NULL,
      subtitle = NULL,
      legend_title = "CD3"
    )
  )
  expect_equal(
    .cell_plot_legend_title(cell_annotation(base_plot, legend_title = "CD3")),
    "CD3"
  )
  expect_equal(.cell_plot_legend_title(base_plot), "marker")
  expect_error(cell_coord_rotate(base_plot, max_degree = 361))
  expect_error(cell_coord_rotate(base_plot, axis = c(0, 0, 0)))
  expect_error(cell_coord_rotate(base_plot, axis = c(0, 1)))
  expect_error(cell_coord_rotate(base_plot, origin = "center"))

  expect_equal(
    cell_illuminate(
      base_plot,
      light_direction = c(0, 3, 4),
      lock_light = TRUE
    )$illuminate[c("light_direction", "lock_light")],
    list(light_direction = c(0, 0.6, 0.8), lock_light = TRUE)
  )
  expect_equal(
    cell_coord_rotate(
      base_plot,
      axis = c(0, 2, 0),
      max_degree = -90,
      origin = "centroid"
    )$coord,
    list(
      type = "rotate",
      axis = c(0, 1, 0),
      max_degree = -90,
      origin = "centroid"
    )
  )

  size_error <- rlang::catch_cnd(
    cell_node_scale_size(base_plot, sizes = c(6, 2))
  )
  expect_equal(size_error$call[[1]], quote(cell_node_scale_size))

  missing_size <- rlang::catch_cnd(
    cell_node_scale_size(cell_plot(plot_data))
  )
  expect_equal(missing_size$call[[1]], quote(cell_node_scale_size))

  missing_z <- rlang::catch_cnd(cell_plot(dplyr::select(plot_data, -z)))
  expect_equal(missing_z$call[[1]], quote(cell_plot))

  expect_error(
    cell_node_scale_size(cell_plot(plot_data), limits = c(1, 5)),
    "requires a column mapped to"
  )
})

test_that("Constant color, size, and alpha work as expected", {
  plot_data <- tibble::tibble(
    x = c(0, 1),
    y = c(1, 0),
    z = c(-1, 1),
    marker = c(0, 2),
    red = c(3, 4)
  )

  constant_plot <- cell_plot(
    plot_data,
    color = "#6699CC",
    size = 3,
    alpha = 0.5
  )
  expect_equal(
    constant_plot$constant,
    list(color = "#6699CC", size = 3, alpha = 0.5)
  )
  expect_equal(
    constant_plot$mapping[c("color", "size", "alpha")],
    list(color = NULL, size = NULL, alpha = NULL)
  )
  expect_equal(
    as.character(summary(constant_plot)),
    c(
      "<cell_plot>",
      "Rows: 2",
      "Mappings: x, y, z, depth, arrange",
      "Constants: color = #6699CC, size = 3, alpha = 0.5",
      "Modifiers: none"
    )
  )

  # Column names win over color names, and bare names are always columns
  expect_equal(cell_plot(plot_data, color = "red")$mapping$color, "red")
  expect_equal(cell_plot(plot_data, color = red)$mapping$color, "red")
  expect_null(cell_plot(plot_data, color = "red")$constant$color)
  expect_equal(
    cell_plot(dplyr::select(plot_data, -red), color = "red")$constant$color,
    "red"
  )
  expect_equal(cell_plot(plot_data, size = 1L)$constant$size, 1)
  expect_equal(cell_plot(plot_data, alpha = 1 / 2)$constant$alpha, 0.5)

  expect_error(
    cell_plot(plot_data, color = "not-a-color"),
    "neither a column"
  )
  expect_error(cell_plot(plot_data, color = "#FF000080"), "fully opaque")
  expect_error(cell_plot(plot_data, color = 3), "one color")
  expect_error(cell_plot(plot_data, color = c("red", "blue")), "one color")
  expect_error(cell_plot(plot_data, size = c(2, 6)), "one number")
  expect_error(cell_plot(plot_data, size = -1), "size")
  expect_error(cell_plot(plot_data, size = Inf), "finite")
  expect_error(cell_plot(plot_data, alpha = 2), "alpha")
  expect_error(cell_plot(plot_data, alpha = TRUE), "one number")

  expect_error(
    cell_plot(plot_data, color = "#6699CC") |>
      cell_node_scale_color(colors = c("blue", "red")),
    "no scale"
  )
  expect_error(
    cell_plot(plot_data, size = 3) |> cell_node_scale_size(),
    "no scale"
  )
  expect_error(
    cell_plot(plot_data, alpha = 0.5) |> cell_node_scale_alpha(),
    "no scale"
  )
})

test_that("Constant node size works as expected", {
  plot_data <- tibble::tibble(
    x = c(0, 1),
    y = c(1, 0),
    z = c(-1, 1)
  )

  expect_equal(
    (cell_plot(plot_data) |> build_cell_plot())$size$resolved,
    1
  )
  expect_equal(
    (cell_plot(plot_data) |> build_cell_plot())$alpha$resolved,
    1
  )
  expect_equal(
    (cell_plot(plot_data) |> build_cell_plot())$color$resolved,
    "#E5E5E5"
  )

  constant_size <- cell_plot(plot_data, size = 2) |>
    build_cell_plot()
  expect_equal(
    constant_size$size,
    list(
      sizes = 2,
      limits = NULL,
      values = NULL,
      resolved = 2
    )
  )

  expect_error(cell_node_scale_size(cell_plot(plot_data), sizes = c(2, 6)))
  expect_error(cell_node_scale_size(cell_plot(plot_data, size = 2)))
  expect_error(cell_node_scale_size(cell_plot(plot_data), sizes = -1))
})

test_that("Categorical size and alpha accept named per-level values", {
  plot_data <- tibble::tibble(
    x = c(0, 1, 2, 3),
    y = c(1, 0, 1, 0),
    z = c(-1, 1, 0, 2),
    group = c("a", "b", "c", NA_character_),
    group_fct = factor(c("a", "b", "c", "a"), levels = c("a", "b", "c"))
  )

  named_sizes <- c(a = 1, b = 2, c = 3)
  sized <- cell_plot(plot_data, size = group) |>
    cell_node_scale_size(sizes = named_sizes)
  expect_equal(sized$size, list(sizes = named_sizes, limits = NULL))

  expect_warning(
    built <- build_cell_plot(sized),
    "Removed 1 row containing missing values in size and/or alpha."
  )
  expect_equal(built$data$group, c("a", "c", "b"))
  expect_equal(unname(built$size$resolved), c(1, 3, 2))
  expect_null(built$size$limits)
  expect_equal(built$size$sizes, named_sizes)

  reversed <- cell_plot(plot_data, size = group_fct, arrange = NULL) |>
    cell_node_scale_size(sizes = c(c = 1, a = 6, b = 2, extra = 9)) |>
    build_cell_plot()
  expect_equal(unname(reversed$size$resolved), c(6, 2, 1, 6))

  expect_warning(
    alpha_built <- cell_plot(plot_data, alpha = group, arrange = NULL) |>
      cell_node_scale_alpha(alphas = c(a = 0.1, b = 0.5, c = 1)) |>
      build_cell_plot(),
    "Removed 1 row containing missing values in size and/or alpha."
  )
  expect_equal(unname(alpha_built$alpha$resolved), c(0.1, 0.5, 1))

  expect_warning(
    ranged <- cell_plot(plot_data, size = group, arrange = NULL) |>
      cell_node_scale_size(sizes = c(2, 6)) |>
      build_cell_plot(),
    "Removed 1 row containing missing values in size and/or alpha."
  )
  expect_equal(unname(ranged$size$resolved), c(2, 4, 6))

  expect_error(
    cell_plot(plot_data, size = group) |>
      cell_node_scale_size(sizes = c(a = 1, b = 2)),
    "names\\(sizes\\)"
  )
  expect_error(
    cell_plot(plot_data, size = group) |>
      cell_node_scale_size(sizes = c(a = 1, 2, c = 3)),
    "name every element"
  )
  expect_error(
    cell_plot(plot_data, size = group) |>
      cell_node_scale_size(
        sizes = c(a = 1, b = 2, c = 3),
        limits = c(1, 2)
      ),
    "cannot be combined with"
  )
  expect_error(
    cell_plot(plot_data, size = z) |>
      cell_node_scale_size(sizes = c(a = 1, b = 2, c = 3))
  )
  expect_error(
    cell_plot(plot_data, size = group) |>
      cell_node_scale_size(sizes = c(1, 2, 3))
  )
  expect_error(
    cell_plot(plot_data, alpha = group) |>
      cell_node_scale_alpha(alphas = c(a = 0.2, b = 1.5, c = 0.4))
  )

  exact_names <- cell_plot(
    dplyr::mutate(plot_data, group = c("a", "alpha", "a", "alpha")),
    size = group,
    arrange = NULL
  ) |>
    cell_node_scale_size(sizes = c(a = 1, alpha = 5)) |>
    build_cell_plot()
  expect_equal(unname(exact_names$size$resolved), c(1, 5, 1, 5))

  unused_level <- cell_plot(
    dplyr::mutate(
      plot_data,
      group = factor(c("a", "b", "a", "b"), levels = c("a", "b", "c"))
    ),
    size = group,
    arrange = NULL
  ) |>
    cell_node_scale_size(sizes = c(a = 1, b = 4))
  expect_equal(
    unused_level$size$sizes,
    c(a = 1, b = 4)
  )
})

test_that("Depth mapping defaults and collisions work as expected", {
  plot_data <- tibble::tibble(
    apa = c(0, 1),
    x = c(2, 3),
    y = c(1, 0),
    z = c(-1, 1)
  )

  default_depth <- cell_plot(plot_data, x = apa)
  expect_equal(default_depth$mapping$depth, "z")
  expect_null(default_depth$depth)

  disabled_depth <- cell_plot(plot_data, x = apa, depth = NULL)
  expect_null(disabled_depth$mapping$depth)

  expect_no_error(cell_plot(plot_data, x = apa, depth = x))
  expect_error(cell_plot(plot_data, x = apa, depth = apa))
  expect_error(cell_plot(plot_data, depth = y))
  expect_error(cell_plot(plot_data, size = z))
  expect_equal(
    cell_plot(plot_data, size = z, depth = NULL)$mapping$size,
    "z"
  )
})

test_that("Palette colors with embedded alpha are rejected", {
  plot_data <- tibble::tibble(x = 0, y = 0, z = 0, marker = 1, cell = "a")
  base_plot <- cell_plot(plot_data, color = marker)
  categorical_plot <- cell_plot(plot_data, color = cell)

  expect_error(
    cell_node_scale_color(base_plot, colors = "#FF000080"),
    "fully opaque"
  )
  expect_error(
    cell_node_scale_color(
      base_plot,
      colors = c("blue", "red"),
      na_color = "#00000080"
    ),
    "fully opaque"
  )
  expect_error(
    cell_node_scale_color(categorical_plot, colors = c(a = "#00FF0080")),
    "fully opaque"
  )
  expect_no_error(
    cell_node_scale_color(base_plot, colors = "#FF0000FF")
  )
  expect_no_error(
    cell_node_scale_color(base_plot, colors = c("black", "white"), limits = 0:1)
  )
})

test_that("Color scale training follows the modifier contract", {
  plot_data <- tibble::tibble(
    x = c(0, 1, 0, 1),
    y = c(0, 0, 1, 1),
    z = c(0, 0, 0, 0),
    marker = c(0, 1, 0, 10),
    cell = c("a", "a", "b", "b")
  )

  trained <- cell_plot(plot_data, color = marker) |>
    cell_grid(cols = cell) |>
    cell_node_scale_color(colors = c("black", "white")) |>
    build_cell_plot()
  expect_equal(trained$color$limits, c(0, 10))
  expect_equal(
    trained$color$resolved,
    .cell_colors_to_hex(
      scales::gradient_n_pal(c("black", "white"))(c(0, 0.1, 0, 1))
    )
  )

  diverging_data <- tibble::tibble(x = 0:2, y = 0:2, z = 0:2, marker = c(-1, 0, 2))
  observed <- cell_plot(diverging_data, color = marker) |>
    cell_node_scale_color(colors = c("blue", "white", "red")) |>
    build_cell_plot()
  expect_equal(observed$color$limits, c(-1, 2))
  expect_equal(
    observed$color$resolved,
    .cell_colors_to_hex(
      scales::gradient_n_pal(c("blue", "white", "red"))(c(0, 1 / 3, 1))
    )
  )

  centered <- cell_plot(diverging_data, color = marker) |>
    cell_node_scale_color(
      colors = c("blue", "white", "red"),
      limits = c(-2, 2)
    ) |>
    build_cell_plot()
  expect_equal(centered$color$limits, c(-2, 2))
  expect_equal(
    centered$color$resolved,
    .cell_colors_to_hex(
      scales::gradient_n_pal(c("blue", "white", "red"))(c(0.25, 0.5, 1))
    )
  )

  categorical_numeric <- cell_plot(
    tibble::tibble(x = 0:1, y = 0:1, z = 0:1, marker = c(2, 1)),
    color = marker
  ) |>
    cell_node_scale_color(
      colors = c("red", "blue"),
      type = "categorical"
    ) |>
    build_cell_plot()
  expect_equal(categorical_numeric$color$type, "categorical")
  expect_equal(categorical_numeric$color$limits, c("1", "2"))
  expect_equal(categorical_numeric$color$resolved, c("#0000FF", "#FF0000"))

  expect_error(
    cell_plot(tibble::tibble(x = 0, y = 0, z = 0, cell = "a"), color = cell) |>
      cell_node_scale_color(colors = "red", type = "continuous")
  )
})

test_that("Relative node sizes convert to backend units and keep depth scaling", {
  expect_equal(.cell_relative_size_to_ggplot(c(2, 6)), c(2, 6))
  expect_equal(.cell_relative_size_to_pixels(c(2, 6)), c(2, 6) * 96 / 25.4)
  expect_equal(
    .cell_relative_size_to_cex(c(2, 6)),
    c(2, 6) / (12 * 25.4 / 72)
  )

  z <- c(0, -10)
  apparent <- .cell_depth_apparent_size(z, base_size = 2, focal_distance = 10)
  expect_equal(apparent, 2 * c(1, 0.5) / (10 / 15))
  expect_equal(apparent[2] / apparent[1], 0.5)
  expect_gt(
    diff(.cell_depth_apparent_size(c(-1, 1), 2, focal_distance = 1)),
    diff(.cell_depth_apparent_size(c(-1, 1), 2, focal_distance = 100))
  )
  expect_equal(
    .cell_depth_apparent_size(c(-1, 0, 1), base_size = 2)[2],
    2
  )
  expect_equal(.cell_depth_apparent_size(c(1, 1), base_size = 2), c(2, 2))
  expect_equal(
    .cell_relative_size_to_ggplot(apparent),
    .cell_relative_size_to_cex(apparent) * (12 * 25.4 / 72)
  )
})

test_that("ggplot depth sizing preserves mean size and ignores mapped size", {
  plot_data <- tibble::tibble(
    x = c(0, 1, 2),
    y = c(0, 0, 0),
    z = c(-1, 0, 1)
  )

  depth_plot <- cell_plot(plot_data, size = 2) |>
    cell_node_depth(focal_distance = 5)
  expect_equal(depth_plot$depth, list(focal_distance = 5))
  depth_built <- build_cell_plot(depth_plot)
  expect_equal(depth_built$size$resolved, 2)
  expect_equal(
    ggplot2::ggplot_build(.render_cell_plot_ggplot(depth_built))$data[[1]]$size,
    .cell_depth_apparent_size(
      c(-1, 0, 1),
      base_size = 2,
      focal_distance = 5
    )
  )

  mapped_built <- cell_plot(plot_data, size = x) |>
    cell_node_scale_size(sizes = c(2, 6)) |>
    build_cell_plot()
  expect_equal(
    ggplot2::ggplot_build(.render_cell_plot_ggplot(mapped_built))$data[[1]]$size,
    mapped_built$size$resolved
  )

  disabled_depth_built <- cell_plot(plot_data, size = 2, depth = NULL) |>
    build_cell_plot()
  expect_equal(
    ggplot2::ggplot_build(
      .render_cell_plot_ggplot(disabled_depth_built)
    )$data[[1]]$size,
    rep(2, 3)
  )

  expect_error(
    cell_plot(plot_data, depth = NULL) |>
      cell_node_depth()
  )
  expect_error(cell_plot(plot_data) |> cell_node_depth(focal_distance = 0))
  expect_error(cell_plot(plot_data) |> cell_node_depth(focal_distance = Inf))
})

test_that("Depth sizing and node size scaling cannot both be used", {
  plot_data <- tibble::tibble(
    x = c(0, 1, 2),
    y = c(0, 0, 0),
    z = c(-1, 0, 1),
    abundance = c(1, 3, 5)
  )

  expect_error(
    cell_plot(plot_data, size = abundance) |>
      cell_node_scale_size(sizes = c(2, 6)) |>
      cell_node_depth(),
    "cannot both be used"
  )
  expect_error(
    cell_plot(plot_data, size = abundance) |>
      cell_node_depth() |>
      cell_node_scale_size(sizes = c(2, 6)),
    "cannot both be used"
  )
  expect_no_error(
    cell_plot(plot_data, size = 3) |>
      cell_node_depth()
  )
  expect_no_error(
    cell_plot(plot_data, size = abundance) |>
      cell_node_scale_size(sizes = c(2, 6))
  )
})
