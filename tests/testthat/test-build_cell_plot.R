test_that("Cell plot building and default rendering work as expected", {
  plot_data <- tibble::tibble(
    x = c(0, 1),
    y = c(1, 0),
    z = c(-1, 1),
    marker = c(0, 2),
    abundance = c(1, 5),
    confidence = c(0, 1),
    cell = c("a", "b")
  )

  plot_recipe <- cell_plot(
    plot_data,
    color = marker,
    size = abundance,
    alpha = confidence
  ) |>
    cell_grid(cols = cell) |>
    cell_node_scale_color(
      colors = c("black", "white"),
      limits = c(0, 2)
    ) |>
    cell_node_scale_size(sizes = c(2, 6), limits = c(1, 5)) |>
    cell_node_scale_alpha(alphas = c(0.2, 1), limits = c(0, 1)) |>
    cell_annotation(title = "Cell layouts")

  built <- build_cell_plot(plot_recipe)

  expect_s3_class(built, "cell_plot_built")
  expect_equal(names(built), names(plot_recipe))
  expect_equal(
    unclass(built),
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
      grid = list(rows = NULL, cols = "cell"),
      color = list(
        colors = c("black", "white"),
        limits = c(0, 2),
        na_color = "grey50",
        type = "continuous",
        values = c(0, 2),
        resolved = c("#000000", "#FFFFFF"),
        illuminated = NULL
      ),
      size = list(
        sizes = c(2, 6),
        limits = c(1, 5),
        values = c(1, 5),
        resolved = c(2, 6)
      ),
      depth = NULL,
      alpha = list(
        alphas = c(0.2, 1),
        limits = c(0, 1),
        values = c(0, 1),
        resolved = c(0.2, 1)
      ),
      theme = list(
        background_color = "white",
        strip_background_color = "#D9D9D9",
        text_color = "black",
        text_size = 11
      ),
      annotation = list(
        title = "Cell layouts",
        subtitle = NULL,
        legend_title = NULL
      ),
      illuminate = NULL,
      coord = NULL
    )
  )

  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  rendered <- print(plot_recipe)
  grDevices::dev.off()
  unlink(plot_file)
  expect_s3_class(rendered, "ggplot")
  expect_s3_class(rendered$facet, "FacetGrid")
  expect_equal(rendered$labels$title, "Cell layouts")
  expect_equal(rendered$scales$get_scales("colour")$name, "marker")
  custom_legend <- cell_plot(plot_data, color = marker) |>
    cell_annotation(legend_title = "CD3") |>
    build_cell_plot()
  expect_equal(
    .render_cell_plot_ggplot(custom_legend)$scales$get_scales("colour")$name,
    "CD3"
  )
  expect_equal(ggplot2::ggplot_build(rendered)$data[[1]]$size, c(2, 6))
  expect_equal(
    unclass(ggplot2::ggplot_build(rendered)$data[[1]]$colour),
    built$color$resolved
  )
  facet_rendered <- cell_plot(
    dplyr::mutate(plot_data, row = cell)
  ) |>
    cell_grid(rows = row, cols = cell) |>
    build_cell_plot() |>
    .render_cell_plot_ggplot()
  expect_equal(
    list(
      switch = facet_rendered$facet$params$switch,
      strip_margin = facet_rendered$theme$strip.text$margin,
      row_strip_angle = facet_rendered$theme$strip.text.y.left$angle
    ),
    list(
      switch = "y",
      strip_margin = ggplot2::margin(4.4, 4.4, 4.4, 4.4),
      row_strip_angle = 90
    )
  )
  grDevices::pdf(NULL)
  facet_grob <- ggplot2::ggplotGrob(facet_rendered)
  top_strip <- facet_grob$layout[
    facet_grob$layout$name == "strip-t-1",
    ,
    drop = FALSE
  ]
  left_strip <- facet_grob$layout[
    facet_grob$layout$name == "strip-l-1",
    ,
    drop = FALSE
  ]
  top_strip_height <- grid::convertHeight(
    sum(facet_grob$heights[top_strip$t:top_strip$b]),
    "pt",
    valueOnly = TRUE
  )
  left_strip_width <- grid::convertWidth(
    sum(facet_grob$widths[left_strip$l:left_strip$r]),
    "pt",
    valueOnly = TRUE
  )
  grDevices::dev.off()
  expect_equal(
    list(
      equal_thickness = isTRUE(all.equal(
        top_strip_height,
        left_strip_width
      )),
      corners_touch = c(
        horizontal = left_strip$r + 1L == top_strip$l,
        vertical = top_strip$b + 1L == left_strip$t
      )
    ),
    list(
      equal_thickness = TRUE,
      corners_touch = c(horizontal = TRUE, vertical = TRUE)
    )
  )
  expect_no_warning(
    .render_cell_plot_ggplot(
      cell_plot(plot_data) |>
        cell_illuminate() |>
        build_cell_plot()
    )
  )

  unordered <- tibble::tibble(
    x = c(0, 1),
    y = c(0, 1),
    z = c(1, -1),
    marker = c(0, 1)
  )
  unordered_built <- cell_plot(unordered, color = marker) |>
    build_cell_plot()
  expect_equal(unordered_built$data$z, c(-1, 1))
  expect_equal(
    ggplot2::ggplot_build(.render_cell_plot_ggplot(unordered_built))$data[[1]]$x,
    c(1, 0)
  )

  categorical_built <- cell_plot(plot_data, color = cell) |>
    cell_node_scale_color(colors = c(a = "red", b = "blue")) |>
    build_cell_plot()
  expect_equal(
    categorical_built$color$resolved,
    c("#FF0000", "#0000FF")
  )
  expect_equal(categorical_built$color$type, "categorical")

  expect_error(cell_plot(plot_data, x = missing_column))
  expect_error(cell_plot(dplyr::mutate(plot_data, x = as.character(x))))
  expect_error(cell_plot(plot_data) |> cell_node_scale_size())
  expect_error(cell_plot(dplyr::select(plot_data, -z)))
  expect_error(cell_plot(plot_data, x = NULL))
  expect_error(cell_plot(plot_data, y = NULL))
  expect_error(cell_plot(plot_data, z = NULL))
  expect_error(
    cell_plot(plot_data, color = cell) |>
      cell_node_scale_color(colors = c("black", "white"), limits = c(0, 1))
  )
})

test_that("Cell illumination resolves colors after scale training", {
  skip_if_not_installed("FNN")

  illumination_data <- tibble::tibble(
    x = rep(c(0, 1, 0), 2),
    y = rep(c(0, 0, 1), 2),
    z = rep(c(-1, 0, 1), 2),
    panel = rep(c("a", "b"), each = 3),
    marker = factor(
      c("signal", "signal", "signal", "signal", NA, "signal"),
      levels = "signal"
    )
  )
  base_recipe <- cell_plot(illumination_data, color = marker) |>
    cell_grid(cols = panel) |>
    cell_node_scale_color(colors = c(signal = "#6699CC"))
  illuminated <- base_recipe |>
    cell_illuminate(
      clamp_quantiles = c(0, 1),
      directional_light_weight = 1,
      volume_shading_weight = 0,
      ambient_occlusion_weight = 0,
      light_direction = c(0, 0, 1)
    ) |>
    build_cell_plot()

  expect_equal(
    illuminated$color,
    list(
      colors = c(signal = "#6699CC"),
      limits = "signal",
      na_color = "grey50",
      type = "categorical",
      values = factor(
        c("signal", "signal", "signal", NA, "signal", "signal"),
        levels = "signal"
      ),
      resolved = c(
        "#6699CC", "#6699CC", "#6699CC", "#7F7F7F", "#6699CC",
        "#6699CC"
      ),
      illuminated = c(
        "#12273D", "#12273D", "#345C85", "#7F7F7F", "#6699CC",
        "#6699CC"
      )
    )
  )

  locked_light <- base_recipe |>
    cell_illuminate(
      clamp_quantiles = c(0, 1),
      directional_light_weight = 1,
      volume_shading_weight = 0,
      ambient_occlusion_weight = 0,
      light_direction = c(1, 0, 0),
      lock_light = TRUE
    ) |>
    build_cell_plot()
  unlocked_light <- base_recipe |>
    cell_illuminate(
      clamp_quantiles = c(0, 1),
      directional_light_weight = 1,
      volume_shading_weight = 0,
      ambient_occlusion_weight = 0,
      light_direction = c(1, 0, 0)
    ) |>
    build_cell_plot()
  expected_light_colors <- list(
    resolved = c(
      "#6699CC", "#6699CC", "#6699CC", "#7F7F7F", "#6699CC",
      "#6699CC"
    ),
    illuminated = c(
      "#12273D", "#12273D", "#6699CC", "#7F7F7F", "#12273D",
      "#12273D"
    )
  )
  expect_equal(
    list(
      locked = locked_light$color[c("resolved", "illuminated")],
      unlocked = unlocked_light$color[c("resolved", "illuminated")]
    ),
    list(
      locked = expected_light_colors,
      unlocked = expected_light_colors
    )
  )

  constant_color <- cell_plot(illumination_data) |>
    cell_grid(cols = panel) |>
    cell_illuminate(
      clamp_quantiles = c(0, 1),
      directional_light_weight = 1,
      volume_shading_weight = 0,
      ambient_occlusion_weight = 0,
      light_direction = c(0, 0, 1)
    ) |>
    build_cell_plot()
  expect_equal(
    constant_color$color,
    list(
      colors = "#E5E5E5",
      limits = NULL,
      na_color = NULL,
      type = NULL,
      values = NULL,
      resolved = "#E5E5E5",
      illuminated = c(
        "#454545", "#454545", "#959595", "#959595", "#E5E5E5",
        "#E5E5E5"
      )
    )
  )
  constant_rendered <- .render_cell_plot_ggplot(constant_color)
  expect_equal(
    list(
      data_colors = as.character(constant_rendered$data$.cell_color),
      rendered_colors = unclass(
        ggplot2::ggplot_build(constant_rendered)$data[[1]]$colour
      )
    ),
    list(
      data_colors = constant_color$color$illuminated,
      rendered_colors = constant_color$color$illuminated
    )
  )

  shadow_color <- base_recipe |>
    cell_illuminate(
      clamp_quantiles = c(0, 1),
      directional_light_weight = 1,
      volume_shading_weight = 0,
      ambient_occlusion_weight = 0,
      shadow_colors = c("black", "blue"),
      light_direction = c(0, 0, 1)
    ) |>
    build_cell_plot()
  expect_equal(
    shadow_color$color$illuminated,
    c(
      "#1E2D3D", "#1E2D3D", "#4E69AE", "#7F7F7F", "#6699CC",
      "#6699CC"
    )
  )

  singleton_panel <- tibble::tibble(
    x = c(0, 1, 0, 0),
    y = c(0, 0, 1, 0),
    z = c(-1, 0, 1, 0),
    panel = c("a", "a", "a", "b"),
    marker = factor(rep("signal", 4))
  ) |>
    cell_plot(color = marker) |>
    cell_grid(cols = panel) |>
    cell_node_scale_color(colors = c(signal = "#6699CC")) |>
    cell_illuminate(
      clamp_quantiles = c(0, 1),
      directional_light_weight = 1,
      volume_shading_weight = 0,
      ambient_occlusion_weight = 0
    ) |>
    build_cell_plot()
  expect_equal(
    singleton_panel$color$illuminated[singleton_panel$data$panel == "b"],
    "#6699CC"
  )

  constant_mask <- tibble::tibble(
    x = c(-1, 0, 1),
    y = 0,
    z = 0,
    marker = factor(rep("signal", 3))
  ) |>
    cell_plot(color = marker) |>
    cell_node_scale_color(colors = c(signal = "#6699CC")) |>
    cell_illuminate(
      clamp_quantiles = c(0, 1),
      directional_light_weight = 1,
      volume_shading_weight = 0,
      ambient_occlusion_weight = 0,
      light_direction = c(0, 0, 1)
    ) |>
    build_cell_plot()
  expect_equal(
    constant_mask$color$illuminated,
    rep("#6699CC", 3)
  )

  expect_error(
    illumination_data |>
      dplyr::select(-z) |>
      cell_plot(color = marker)
  )
  expect_error(
    base_recipe |>
      cell_illuminate(
        directional_light_weight = 0,
        volume_shading_weight = 0,
        ambient_occlusion_weight = 0
      ) |>
      build_cell_plot()
  )
})

test_that("Constant color, size, and alpha reach every rendered point", {
  plot_data <- tibble::tibble(
    x = c(0, 1, 2),
    y = c(0, 1, 2),
    z = c(-1, 0, 1),
    marker = c(0, 1, 2)
  )

  built <- cell_plot(
    plot_data,
    color = "red",
    size = 3,
    alpha = 0.5,
    depth = NULL
  ) |>
    build_cell_plot()

  expect_equal(
    built$constant,
    list(color = "red", size = 3, alpha = 0.5)
  )
  expect_equal(
    list(color = built$color, size = built$size, alpha = built$alpha),
    list(
      color = list(
        colors = "#FF0000",
        limits = NULL,
        na_color = NULL,
        type = NULL,
        values = NULL,
        resolved = "#FF0000",
        illuminated = NULL
      ),
      size = list(sizes = 3, limits = NULL, values = NULL, resolved = 3),
      alpha = list(alphas = 0.5, limits = NULL, values = NULL, resolved = 0.5)
    )
  )

  rendered <- ggplot2::ggplot_build(.render_cell_plot_ggplot(built))$data[[1]]
  expect_equal(
    list(
      color = unclass(rendered$colour),
      size = rendered$size,
      alpha = rendered$alpha
    ),
    list(
      color = rep("#FF0000", 3),
      size = rep(3, 3),
      alpha = rep(0.5, 3)
    )
  )
  expect_null(
    .render_cell_plot_ggplot(built)$scales$get_scales("colour")
  )

  depth_built <- cell_plot(plot_data, size = 3) |>
    build_cell_plot()
  expect_equal(depth_built$size$resolved, 3)
  expect_equal(
    ggplot2::ggplot_build(.render_cell_plot_ggplot(depth_built))$data[[1]]$size,
    .cell_depth_apparent_size(c(-1, 0, 1), base_size = 3)
  )

  constant_alpha_colors <- cell_plot(plot_data, color = "red", alpha = 0.5) |>
    build_cell_plot()
  expect_equal(
    scales::alpha(
      rep_len(.cell_plot_rendered_colors(constant_alpha_colors), 3),
      alpha = rep_len(constant_alpha_colors$alpha$resolved, 3)
    ),
    rep("#FF000080", 3)
  )
})

test_that("Color palettes are stable and match ggplot interpolation", {
  clipped_color <- cell_plot(
    tibble::tibble(x = 0:2, y = 0:2, z = 0:2, marker = c(-1, 1, 3)),
    color = marker
  ) |>
    cell_node_scale_color(colors = c("black", "white"), limits = c(0, 2)) |>
    build_cell_plot()
  expect_equal(
    clipped_color$color$resolved,
    .cell_colors_to_hex(
      scales::gradient_n_pal(c("black", "white"))(c(0, 0.5, 1))
    )
  )
  rendered_limits <- .render_cell_plot_ggplot(clipped_color)
  color_scale <- rendered_limits$scales$get_scales("colour")
  expect_equal(color_scale$oob(c(-1, 3), range = c(0, 2)), c(0, 2))
  rendered_colors <- ggplot2::ggplot_build(rendered_limits)$data[[1]]$colour
  expect_equal(
    unclass(substr(rendered_colors, 1, 7)),
    clipped_color$color$resolved
  )

  palette_order <- tibble::tibble(
    x = c(0, 1, 2),
    y = c(0, 1, 2),
    z = c(3, 1, 2),
    group = c("b", "a", "c")
  )
  colors_by_z <- cell_plot(palette_order, color = group, arrange = z) |>
    cell_node_scale_color(colors = c("red", "green", "blue")) |>
    build_cell_plot()
  colors_by_x <- cell_plot(palette_order, color = group, arrange = x) |>
    cell_node_scale_color(colors = c("red", "green", "blue")) |>
    build_cell_plot()
  resolved_by_group <- function(built) {
    resolved <- built$color$resolved
    names(resolved) <- built$data$group
    return(unname(resolved[c("a", "b", "c")]))
  }
  expect_equal(resolved_by_group(colors_by_z), c("#FF0000", "#00FF00", "#0000FF"))
  expect_equal(resolved_by_group(colors_by_z), resolved_by_group(colors_by_x))
  expect_equal(colors_by_z$color$limits, c("a", "b", "c"))
  expect_equal(.cell_colors_to_hex("red"), "#FF0000")
})

test_that("Color, size, and alpha mappings accept character, factor, numeric, and integer columns", {
  plot_data <- tibble::tibble(
    x = c(0, 1),
    y = c(1, 0),
    z = c(-1, 1),
    color_chr = c("a", "b"),
    color_fct = factor(c("a", "b"), levels = c("a", "b")),
    color_dbl = c(0, 1),
    color_int = c(0L, 1L),
    size_chr = c("a", "b"),
    size_fct = factor(c("a", "b"), levels = c("a", "b")),
    size_dbl = c(1, 5),
    size_int = c(1L, 5L),
    alpha_chr = c("a", "b"),
    alpha_fct = factor(c("a", "b"), levels = c("a", "b")),
    alpha_dbl = c(0, 1),
    alpha_int = c(0L, 1L)
  )

  categorical_color <- list(
    colors = c(a = "red", b = "blue"),
    limits = c("a", "b"),
    na_color = "grey50",
    type = "categorical",
    values = c("a", "b"),
    resolved = c("#FF0000", "#0000FF"),
    illuminated = NULL
  )
  continuous_color <- list(
    colors = c("black", "white"),
    limits = c(0, 1),
    na_color = "grey50",
    type = "continuous",
    values = c(0, 1),
    resolved = c("#000000", "#FFFFFF"),
    illuminated = NULL
  )
  discrete_size <- list(
    sizes = c(2, 6),
    limits = c(1, 2),
    values = c("a", "b"),
    resolved = c(2, 6)
  )
  continuous_size <- list(
    sizes = c(2, 6),
    limits = c(1, 5),
    values = c(1, 5),
    resolved = c(2, 6)
  )
  discrete_alpha <- list(
    alphas = c(0.2, 1),
    limits = c(1, 2),
    values = c("a", "b"),
    resolved = c(0.2, 1)
  )
  continuous_alpha <- list(
    alphas = c(0.2, 1),
    limits = c(0, 1),
    values = c(0, 1),
    resolved = c(0.2, 1)
  )

  character_built <- cell_plot(
    plot_data,
    color = color_chr,
    size = size_chr,
    alpha = alpha_chr
  ) |>
    cell_node_scale_color(colors = c(a = "red", b = "blue")) |>
    cell_node_scale_size(sizes = c(2, 6)) |>
    cell_node_scale_alpha(alphas = c(0.2, 1)) |>
    build_cell_plot()
  expect_equal(
    list(
      color = character_built$color,
      size = character_built$size,
      alpha = character_built$alpha
    ),
    list(
      color = categorical_color,
      size = discrete_size,
      alpha = discrete_alpha
    )
  )

  factor_built <- cell_plot(
    plot_data,
    color = color_fct,
    size = size_fct,
    alpha = alpha_fct
  ) |>
    cell_node_scale_color(colors = c(a = "red", b = "blue")) |>
    cell_node_scale_size(sizes = c(2, 6)) |>
    cell_node_scale_alpha(alphas = c(0.2, 1)) |>
    build_cell_plot()
  categorical_factor_color <- categorical_color
  categorical_factor_color$values <- plot_data$color_fct
  discrete_factor_size <- discrete_size
  discrete_factor_size$values <- plot_data$size_fct
  discrete_factor_alpha <- discrete_alpha
  discrete_factor_alpha$values <- plot_data$alpha_fct
  expect_equal(
    list(
      color = factor_built$color,
      size = factor_built$size,
      alpha = factor_built$alpha
    ),
    list(
      color = categorical_factor_color,
      size = discrete_factor_size,
      alpha = discrete_factor_alpha
    )
  )

  numeric_built <- cell_plot(
    plot_data,
    color = color_dbl,
    size = size_dbl,
    alpha = alpha_dbl
  ) |>
    cell_node_scale_color(colors = c("black", "white"), limits = c(0, 1)) |>
    cell_node_scale_size(sizes = c(2, 6), limits = c(1, 5)) |>
    cell_node_scale_alpha(alphas = c(0.2, 1), limits = c(0, 1)) |>
    build_cell_plot()
  expect_equal(
    list(
      color = numeric_built$color,
      size = numeric_built$size,
      alpha = numeric_built$alpha
    ),
    list(
      color = continuous_color,
      size = continuous_size,
      alpha = continuous_alpha
    )
  )

  integer_built <- cell_plot(
    plot_data,
    color = color_int,
    size = size_int,
    alpha = alpha_int
  ) |>
    cell_node_scale_color(colors = c("black", "white"), limits = c(0, 1)) |>
    cell_node_scale_size(sizes = c(2, 6), limits = c(1, 5)) |>
    cell_node_scale_alpha(alphas = c(0.2, 1), limits = c(0, 1)) |>
    build_cell_plot()
  continuous_integer_color <- continuous_color
  continuous_integer_color$values <- plot_data$color_int
  continuous_integer_size <- continuous_size
  continuous_integer_size$values <- plot_data$size_int
  continuous_integer_alpha <- continuous_alpha
  continuous_integer_alpha$values <- plot_data$alpha_int
  expect_equal(
    list(
      color = integer_built$color,
      size = integer_built$size,
      alpha = integer_built$alpha
    ),
    list(
      color = continuous_integer_color,
      size = continuous_integer_size,
      alpha = continuous_integer_alpha
    )
  )

  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  expect_s3_class(print(cell_plot(plot_data, color = color_chr, size = size_chr, alpha = alpha_chr)), "ggplot")
  expect_s3_class(print(cell_plot(plot_data, color = color_fct, size = size_fct, alpha = alpha_fct)), "ggplot")
  expect_s3_class(print(cell_plot(plot_data, color = color_dbl, size = size_dbl, alpha = alpha_dbl)), "ggplot")
  expect_s3_class(print(cell_plot(plot_data, color = color_int, size = size_int, alpha = alpha_int)), "ggplot")
  grDevices::dev.off()
  unlink(plot_file)
})

test_that("build_cell_plot drops rows with NA in size or alpha with a warning", {
  test_df <- tibble::tibble(
    x = c(1, 2, 3, 4),
    y = c(1, 2, 3, 4),
    z = c(1, 2, 3, 4),
    s = c(1, NA, 2, 3),
    a = c(0.1, 0.2, NA, 0.5),
    c = c("red", "blue", "green", "yellow")
  )

  # Size NA only
  expect_warning(
    built_size <- cell_plot(test_df, size = s) |> build_cell_plot(),
    "Removed 1 row containing missing values in size and/or alpha."
  )
  expect_equal(nrow(built_size$data), 3)
  expect_equal(length(built_size$size$resolved), 3)
  expect_equal(length(built_size$size$values), 3)
  expect_equal(built_size$data$s, c(1, 2, 3))

  # Alpha NA only
  expect_warning(
    built_alpha <- cell_plot(test_df, alpha = a) |> build_cell_plot(),
    "Removed 1 row containing missing values in size and/or alpha."
  )
  expect_equal(nrow(built_alpha$data), 3)
  expect_equal(length(built_alpha$alpha$resolved), 3)
  expect_equal(length(built_alpha$alpha$values), 3)
  expect_equal(built_alpha$data$a, c(0.1, 0.2, 0.5))

  # Both size and alpha NAs
  expect_warning(
    built_both <- cell_plot(test_df, size = s, alpha = a, color = c) |> build_cell_plot(),
    "Removed 2 rows containing missing values in size and/or alpha."
  )
  expect_equal(nrow(built_both$data), 2)
  expect_equal(length(built_both$size$resolved), 2)
  expect_equal(length(built_both$size$values), 2)
  expect_equal(length(built_both$alpha$resolved), 2)
  expect_equal(length(built_both$alpha$values), 2)
  expect_equal(length(built_both$color$resolved), 2)
  expect_equal(length(built_both$color$values), 2)
  expect_equal(built_both$data$x, c(1, 4))

  # Illumination is computed after the drop, matching data with NAs already removed
  expect_warning(
    built_illum <- cell_plot(test_df, size = s, color = c) |>
      cell_illuminate() |>
      build_cell_plot(),
    "Removed 1 row containing missing values in size and/or alpha."
  )
  kept_illum <- cell_plot(test_df[!is.na(test_df$s), ], size = s, color = c) |>
    cell_illuminate() |>
    build_cell_plot()
  expect_equal(nrow(built_illum$data), 3)
  expect_equal(length(built_illum$color$illuminated), 3)
  expect_equal(built_illum$color$illuminated, kept_illum$color$illuminated)

  color_after_drop <- tibble::tibble(
    x = c(1, 2),
    y = c(1, 2),
    z = c(1, 2),
    s = c(NA_real_, 1),
    c = c(1, NA_real_)
  )
  expect_warning(
    expect_error(
      cell_plot(color_after_drop, size = s, color = c) |> build_cell_plot(),
      "no usable values"
    ),
    "Removed 1 row"
  )

  pair_df <- tibble::tibble(
    x = c(0, 1),
    y = c(0, 1),
    z = c(0, 1),
    s = c(1, NA_real_),
    c = c("red", "blue")
  )
  expect_warning(
    built_pair <- cell_plot(pair_df, size = s, color = c) |>
      cell_illuminate() |>
      build_cell_plot(),
    "Removed 1 row containing missing values in size and/or alpha."
  )
  built_one <- cell_plot(pair_df[1, , drop = FALSE], size = s, color = c) |>
    cell_illuminate() |>
    build_cell_plot()
  expect_equal(nrow(built_pair$data), 1)
  expect_equal(built_pair$color$illuminated, built_one$color$illuminated)

  # Abort when all rows are dropped across size and alpha
  cross_na_df <- tibble::tibble(
    x = c(1, 2),
    y = c(1, 2),
    z = c(1, 2),
    s = c(1, NA),
    a = c(NA, 0.5)
  )
  expect_error(
    cell_plot(cross_na_df, size = s, alpha = a) |> build_cell_plot(),
    "All rows were dropped due to missing values in size and/or alpha."
  )
})
