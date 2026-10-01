rgl_point_summary <- function(object_id, scene) {
  obj <- scene$objects[[as.character(object_id)]]
  rgba <- rgl::rgl.attrib(as.integer(object_id), "colors")
  return(list(
    type = obj$type,
    x = as.numeric(obj$vertices[, "x"]),
    y = as.numeric(obj$vertices[, "y"]),
    z = as.numeric(obj$vertices[, "z"]),
    size = as.numeric(obj$material$size),
    color = unname(as.character(
      if (is.null(obj$material$color)) character() else obj$material$color
    )),
    rgba = data.frame(
      r = as.numeric(rgba[, "r"]),
      g = as.numeric(rgba[, "g"]),
      b = as.numeric(rgba[, "b"]),
      a = as.numeric(rgba[, "a"])
    )
  ))
}

rgl_all_point_ids <- function(scene) {
  return(as.integer(names(Filter(
    function(obj) identical(obj$type, "points"),
    scene$objects
  ))))
}

rgl_subscene_point_ids <- function(subscene, point_ids) {
  found <- intersect(as.integer(subscene$objects), point_ids)
  for (child in subscene$subscenes) {
    found <- c(found, rgl_subscene_point_ids(child, point_ids))
  }
  return(found)
}

rgl_subscene_object_types <- function(subscene, scene) {
  object_ids <- intersect(
    as.integer(subscene$objects),
    as.integer(names(scene$objects))
  )
  return(vapply(
    as.character(object_ids),
    function(id) scene$objects[[id]]$type,
    character(1),
    USE.NAMES = FALSE
  ))
}

rgl_is_data_panel <- function(subscene, scene) {
  modes <- unname(as.character(subscene$par3d$mouseMode))
  return("trackball" %in% modes)
}

rgl_collect_data_panels <- function(subscene, scene) {
  if (rgl_is_data_panel(subscene, scene) && length(subscene$subscenes) == 0) {
    return(list(subscene))
  }
  if (length(subscene$subscenes) == 0) {
    return(list())
  }
  return(unlist(
    lapply(
      subscene$subscenes,
      rgl_collect_data_panels,
      scene = scene
    ),
    recursive = FALSE
  ))
}

rgl_plot_summary <- function(scene) {
  if (!inherits(scene, "rglscene")) {
    scene <- rgl::scene3d()
  }
  root <- scene$rootSubscene
  point_ids <- rgl_all_point_ids(scene)
  panels <- rgl_collect_data_panels(root, scene)
  if (length(panels) == 0) {
    panels <- list(root)
  }

  panel_summaries <- lapply(panels, function(panel) {
    panel_point_ids <- rgl_subscene_point_ids(panel, point_ids)
    return(list(
      id = panel$id,
      viewport = as.numeric(panel$par3d$viewport),
      zoom = panel$par3d$zoom,
      mouseMode = unname(as.character(panel$par3d$mouseMode)),
      listeners = as.integer(panel$par3d$listeners),
      bbox = as.numeric(panel$par3d$bbox),
      points = lapply(panel_point_ids, rgl_point_summary, scene = scene)
    ))
  })

  backgrounds <- Filter(
    function(obj) identical(obj$type, "background"),
    scene$objects
  )
  background_colors <- unname(vapply(
    backgrounds,
    function(obj) obj$material$color[[1]],
    character(1)
  ))

  object_types <- sort(unname(vapply(
    scene$objects,
    function(obj) obj$type,
    character(1)
  )))

  return(list(
    n_panels = length(panel_summaries),
    panels = panel_summaries,
    background_colors = background_colors,
    windowRect = as.numeric(root$par3d$windowRect),
    n_root_subscenes = length(root$subscenes),
    object_types = object_types
  ))
}

test_that("rgl cell plots work as expected", {
  old_options <- options(rgl.useNULL = TRUE)
  on.exit(options(old_options), add = TRUE)
  on.exit(try(rgl::close3d(), silent = TRUE), add = TRUE)

  plot_data <- tibble::tibble(
    x = c(0, 1),
    y = c(1, 0),
    z = c(-1, 1),
    marker = c(0, 2),
    abundance = c(1, 5),
    confidence = c(0.2, 1),
    cell = c("b", "a")
  )

  continuous <- cell_plot(
    plot_data,
    color = marker,
    size = abundance,
    alpha = confidence
  ) |>
    cell_node_scale_color(colors = c("black", "white"), limits = c(0, 2)) |>
    cell_node_scale_size(sizes = c(2, 6), limits = c(1, 5)) |>
    cell_node_scale_alpha(alphas = c(0.2, 1), limits = c(0, 1)) |>
    cell_theme(
      background_color = "navy",
      text_color = "white",
      text_size = 12
    ) |>
    cell_annotation(title = "Cells", subtitle = "demo") |>
    cell_coord_rotate() |>
    cell_illuminate() |>
    cell_plot_rgl()

  expect_equal(continuous, as.integer(rgl::cur3d()))
  continuous_summary <- rgl_plot_summary(continuous)
  expect_equal(continuous_summary$n_panels, 1L)
  expect_equal(continuous_summary$n_root_subscenes, 3L)
  expect_equal(continuous_summary$windowRect, c(100, 100, 1100, 1100))
  expect_equal(
    continuous_summary$panels[[1]]$points,
    list(
      list(
        type = "points",
        x = 0,
        y = 1,
        z = -1,
        size = 2 * 96 / 25.4,
        color = character(),
        rgba = data.frame(r = 0, g = 0, b = 0, a = 0.356862753629684)
      ),
      list(
        type = "points",
        x = 1,
        y = 0,
        z = 1,
        size = 6 * 96 / 25.4,
        color = "#FFFFFF",
        rgba = data.frame(r = 1, g = 1, b = 1, a = 1)
      )
    ),
    tolerance = 1e-6
  )

  rgl::close3d()
  constant <- cell_plot(plot_data) |>
    cell_plot_rgl()
  constant_summary <- rgl_plot_summary(constant)
  expect_equal(
    constant_summary$panels[[1]]$points[[1]][c("color", "size", "rgba")],
    list(
      color = c("#E5E5E5", "#E5E5E5"),
      size = 1 * 96 / 25.4,
      rgba = data.frame(
        r = c(229, 229) / 255,
        g = c(229, 229) / 255,
        b = c(229, 229) / 255,
        a = c(1, 1)
      )
    ),
    tolerance = 1e-6
  )
  expect_equal(constant_summary$background_colors, "#FFFFFF")

  rgl::close3d()
  mixed <- cell_plot(plot_data, color = cell) |>
    cell_node_scale_color(colors = c(a = "red", b = "blue")) |>
    cell_plot_rgl()
  expect_equal(
    rgl_plot_summary(mixed)$panels[[1]]$points[[1]][c("x", "y", "z", "rgba")],
    list(
      x = c(0, 1),
      y = c(1, 0),
      z = c(-1, 1),
      rgba = data.frame(
        r = c(0, 1),
        g = c(0, 0),
        b = c(1, 0),
        a = c(1, 1)
      )
    )
  )

  rgl::close3d()
  categorical <- cell_plot(plot_data, color = cell) |>
    cell_node_scale_color(colors = c(a = "red", b = "blue")) |>
    cell_grid(cols = cell) |>
    cell_plot_rgl()
  categorical_summary <- rgl_plot_summary(categorical)
  expect_equal(categorical_summary$n_panels, 2L)
  expect_equal(
    lapply(categorical_summary$panels, function(panel) {
      list(
        x = panel$points[[1]]$x,
        rgba = panel$points[[1]]$rgba
      )
    }),
    list(
      list(
        x = 1,
        rgba = data.frame(r = 1, g = 0, b = 0, a = 1)
      ),
      list(
        x = 0,
        rgba = data.frame(r = 0, g = 0, b = 1, a = 1)
      )
    )
  )

  expect_error(cell_plot(dplyr::select(plot_data, -z)))
})

test_that("rgl cell plots ignore depth sizing and arrangement", {
  old_options <- options(rgl.useNULL = TRUE)
  on.exit(options(old_options), add = TRUE)
  on.exit(try(rgl::close3d(), silent = TRUE), add = TRUE)

  plot_data <- tibble::tibble(
    x = c(1, 0),
    y = c(1, 0),
    z = c(-1, 1)
  )

  without_depth <- cell_plot(plot_data, size = 2, depth = NULL) |>
    cell_plot_rgl()
  without_depth <- rgl_plot_summary(without_depth)
  rgl::close3d()
  with_depth <- cell_plot(plot_data, size = 2) |>
    cell_node_depth(focal_distance = 0.2) |>
    cell_plot_rgl()
  with_depth <- rgl_plot_summary(with_depth)
  rgl::close3d()
  arranged <- cell_plot(plot_data, arrange = x) |>
    cell_plot_rgl()
  arranged <- rgl_plot_summary(arranged)

  expect_equal(
    without_depth$panels[[1]]$points[[1]]$size,
    2 * 96 / 25.4,
    tolerance = 1e-6
  )
  expect_equal(
    with_depth$panels[[1]]$points[[1]]$size,
    without_depth$panels[[1]]$points[[1]]$size
  )
  expect_equal(
    arranged$panels[[1]]$points[[1]]$x,
    c(1, 0)
  )
})

test_that("rgl cell plot grids work as expected", {
  old_options <- options(rgl.useNULL = TRUE)
  on.exit(options(old_options), add = TRUE)
  on.exit(try(rgl::close3d(), silent = TRUE), add = TRUE)

  na_factor <- tibble::tibble(
    x = 1:3,
    y = 1:3,
    z = 1:3,
    panel = factor(c("b", NA, "a"), levels = c("a", "b", "unused"))
  ) |>
    cell_plot() |>
    cell_grid(rows = panel) |>
    cell_plot_rgl()

  na_summary <- rgl_plot_summary(na_factor)
  expect_equal(na_summary$n_panels, 4L)
  expect_equal(
    lapply(na_summary$panels, function(panel) {
      if (length(panel$points) == 0) {
        return(list(x = numeric(), n_points = 0L))
      }
      return(list(
        x = panel$points[[1]]$x,
        n_points = length(panel$points[[1]]$x)
      ))
    }),
    list(
      list(x = 3, n_points = 1L),
      list(x = 1, n_points = 1L),
      list(x = numeric(), n_points = 0L),
      list(x = 2, n_points = 1L)
    )
  )
  expect_equal(
    unique(lapply(na_summary$panels, function(panel) round(panel$bbox, 6))),
    list(c(0.9, 3.1, 0.9, 3.1, 0.9, 3.1))
  )

  too_many_rows <- tibble::tibble(
    x = 1:11,
    y = 1:11,
    z = 1:11,
    panel = 1:11
  )
  expect_error(
    cell_plot(too_many_rows) |>
      cell_grid(rows = panel) |>
      cell_plot_rgl()
  )
  too_many_cols <- tibble::tibble(
    x = 1:21,
    y = 1:21,
    z = 1:21,
    panel = 1:21
  )
  expect_error(
    cell_plot(too_many_cols) |>
      cell_grid(cols = panel) |>
      cell_plot_rgl()
  )

  rgl::close3d()
  one_level <- tibble::tibble(
    x = c(0, 1),
    y = c(0, 1),
    z = c(0, 1),
    panel = "only"
  ) |>
    cell_plot() |>
    cell_grid(cols = panel) |>
    cell_plot_rgl()
  expect_equal(
    list(
      n_panels = rgl_plot_summary(one_level)$n_panels,
      n_root_subscenes = rgl_plot_summary(one_level)$n_root_subscenes
    ),
    list(n_panels = 1L, n_root_subscenes = 2L)
  )
})

test_that("rgl cell plot facet chrome and mouse sharing work as expected", {
  old_options <- options(rgl.useNULL = TRUE)
  on.exit(options(old_options), add = TRUE)
  on.exit(try(rgl::close3d(), silent = TRUE), add = TRUE)

  plot_data <- tibble::tibble(
    x = c(0, 1, 2),
    y = c(0, 1, 2),
    z = c(0, 1, 2),
    row = factor(c("r2", "r1", "r2"), levels = c("r1", "r2")),
    col = factor(c("c2", "c2", "c1"), levels = c("c1", "c2")),
    marker = c(0, 1, 2)
  )

  faceted <- cell_plot(plot_data, color = marker) |>
    cell_node_scale_color(colors = c("black", "white"), limits = c(0, 2)) |>
    cell_grid(rows = row, cols = col) |>
    cell_annotation(title = "Facet title") |>
    cell_plot_rgl()

  summary <- rgl_plot_summary(faceted)
  expect_equal(summary$n_panels, 4L)
  # Four dedicated facet strips, the strip-corner cell, a title, and a legend
  # sit outside the four data panels, so points cannot cover facet text.
  expect_equal(summary$n_root_subscenes, 11L)
  expect_equal(
    lapply(summary$panels, function(panel) {
      if (length(panel$points) == 0) {
        return(list(n_points = 0L, x = numeric()))
      }
      return(list(
        n_points = length(panel$points[[1]]$x),
        x = panel$points[[1]]$x
      ))
    }),
    list(
      list(n_points = 0L, x = numeric()),
      list(n_points = 1L, x = 1),
      list(n_points = 1L, x = 2),
      list(n_points = 1L, x = 0)
    )
  )

  panel_ids <- vapply(summary$panels, function(panel) panel$id, numeric(1))
  expect_equal(
    lapply(summary$panels, function(panel) sort(panel$listeners)),
    rep(list(sort(as.integer(panel_ids))), 4)
  )
  # layout3d() chrome must use mouseMode="replace"; otherwise disabling
  # title/legend mouse writes through the inherited parent and leaves every
  # data panel with mouseMode all "none" (non-interactive spin/zoom).
  expect_true(all(vapply(
    summary$panels,
    function(panel) {
      all(c("trackball", "zoom") %in% panel$mouseMode) &&
        !all(panel$mouseMode == "none")
    },
    logical(1)
  )))
  # Title, strips, and the legend are orthographic and ignore the pointer.
  # They do not listen for the data panels, so dragging them does not move
  # the plot.
  all_subscenes <- rgl::scene3d()$rootSubscene$subscenes
  chrome_subscenes <- Filter(
    function(subscene) {
      modes <- unname(as.character(subscene$par3d$mouseMode))
      return(!"trackball" %in% modes)
    },
    all_subscenes
  )
  expect_equal(length(chrome_subscenes), 7L)
  expect_true(all(vapply(
    chrome_subscenes,
    function(subscene) {
      modes <- unname(as.character(subscene$par3d$mouseMode))
      return(all(modes == "none") && identical(subscene$par3d$FOV, 0))
    },
    logical(1)
  )))
  expect_equal(
    lapply(chrome_subscenes, function(subscene) {
      as.integer(subscene$par3d$listeners)
    }),
    lapply(chrome_subscenes, function(subscene) as.integer(subscene$id))
  )
  expect_equal(
    unique(lapply(summary$panels, function(panel) panel$zoom)),
    list(1)
  )
  # Data markers stay the unlit points primitive. Chrome contributes text
  # and the continuous legend's quads, with no bgplot3d textures.
  expect_equal(
    c(
      any(summary$object_types == "points"),
      any(summary$object_types == "text"),
      any(summary$object_types == "quads"),
      any(summary$object_types %in% c("spheres", "sprites", "mesh3d"))
    ),
    c(TRUE, TRUE, TRUE, FALSE)
  )
  expect_false(any(vapply(
    rgl::scene3d()$objects,
    function(object) {
      texture <- object$material$texture
      return(is.character(texture) && any(nzchar(texture)))
    },
    logical(1)
  )))

  expect_equal(
    pixelatorR:::.cell_plot_default_theme$strip_background_color,
    "#D9D9D9"
  )
  both_strips <- pixelatorR:::.cell_rgl_facet_layout(
    n_row = 2L,
    n_col = 2L,
    need_col_strips = TRUE,
    need_row_strips = TRUE,
    need_title = TRUE,
    need_legend = TRUE
  )
  expect_equal(
    both_strips,
    list(
      mat = matrix(
        c(
          10L, 7L, 8L, 9L,
          10L, 5L, 1L, 3L,
          10L, 6L, 2L, 4L,
          11L, 11L, 11L, 11L
        ),
        nrow = 4L
      ),
      widths = c(0.107547169811321, 1, 1, 0.28),
      heights = c(0.12, 0.1, 1, 1)
    )
  )
  rectangular_strips <- pixelatorR:::.cell_rgl_facet_layout(
    n_row = 2L,
    n_col = 2L,
    need_col_strips = TRUE,
    need_row_strips = TRUE,
    need_title = TRUE,
    need_legend = TRUE,
    viewport = c(width = 1200, height = 800)
  )
  expect_equal(
    1200 * rectangular_strips$widths[[1]] /
      sum(rectangular_strips$widths),
    800 * rectangular_strips$heights[[2]] /
      sum(rectangular_strips$heights)
  )
  expect_equal(pixelatorR:::.cell_plot_row_strip_angle, 90)
})

test_that("rgl cell plots use builder-baked illumination colors", {
  old_options <- options(rgl.useNULL = TRUE)
  on.exit(options(old_options), add = TRUE)
  on.exit(try(rgl::close3d(), silent = TRUE), add = TRUE)

  set.seed(1)
  plot_data <- tidyr::expand_grid(
    cell_id = c("c1", "c2"),
    marker = c("m1", "m2"),
    i = 1:30
  ) |>
    dplyr::mutate(
      x = cos(.data$i / 30 * 2 * pi) + rnorm(dplyr::n(), sd = 0.05),
      y = sin(.data$i / 30 * 2 * pi) + rnorm(dplyr::n(), sd = 0.05),
      z = cos(.data$i / 15 * pi) + rnorm(dplyr::n(), sd = 0.05),
      val = 1
    )

  recipe <- cell_plot(plot_data, color = val) |>
    cell_grid(cols = "cell_id", rows = "marker") |>
    cell_node_scale_color(colors = c("lightgrey", "lightgray")) |>
    cell_illuminate(
      ambient_occlusion_weight = 0,
      directional_light_weight = 1,
      volume_shading_weight = 0
    ) |>
    cell_annotation(title = "Spectral layout", subtitle = "demo")

  built <- build_cell_plot(recipe)
  rendered <- pixelatorR:::.cell_plot_rendered_colors(built)
  resolved <- built$color$resolved

  # Flat palette: scale colors are constant; illumination must create shading.
  expect_equal(length(unique(resolved)), 1L)
  expect_equal(length(unique(rendered)) > 1L, TRUE)
  expect_equal(rendered, built$color$illuminated)

  scene <- cell_plot_rgl(recipe)
  summary <- rgl_plot_summary(scene)
  scene_rgba <- do.call(
    rbind,
    lapply(summary$panels, function(panel) {
      do.call(
        rbind,
        lapply(panel$points, function(pt) pt$rgba)
      )
    })
  )
  scene_hex <- grDevices::rgb(
    scene_rgba$r,
    scene_rgba$g,
    scene_rgba$b,
    maxColorValue = 1
  )

  expect_equal(sort(unique(scene_hex)), sort(unique(rendered)))
  expect_equal(any(scene_hex != resolved[[1]]), TRUE)
  # Spin/zoom must remain interactive after the chrome mouseMode fix.
  expect_true(all(vapply(
    summary$panels,
    function(panel) all(c("trackball", "zoom") %in% panel$mouseMode),
    logical(1)
  )))

  # Plotly consumes the same rendered colors for the same built recipe.
  expect_equal(
    pixelatorR:::.cell_plot_rendered_colors(built),
    built$color$illuminated
  )
})

test_that("rgl helpers tolerate missing sizes, flat colorbars, and NA facets", {
  old_options <- options(rgl.useNULL = TRUE)
  on.exit(options(old_options), add = TRUE)
  on.exit(try(rgl::close3d(), silent = TRUE), add = TRUE)

  rgl::open3d()
  ids <- pixelatorR:::.cell_rgl_points(
    x = c(0, 1, 2),
    y = c(0, 1, 2),
    z = c(0, 1, 2),
    color = c("#000000", "#111111", "#222222"),
    size = c(3, NA_real_, 3),
    alpha = c(1, 1, 1)
  )
  expect_true(length(ids) >= 1L)
  rgl::close3d()

  rgl::open3d()
  grouped_ids <- pixelatorR:::.cell_rgl_points(
    x = 1:100,
    y = 1:100,
    z = 1:100,
    color = rep("#000000", 100),
    size = 1:100,
    alpha = rep(1, 100)
  )
  expect_equal(
    list(
      n_objects = length(grouped_ids),
      n_points = sum(vapply(
        grouped_ids,
        function(id) nrow(rgl::scene3d()$objects[[as.character(id)]]$vertices),
        integer(1)
      ))
    ),
    list(n_objects = 20L, n_points = 100L)
  )
  rgl::close3d()

  rgl::open3d()
  pixel_sizes <- pixelatorR:::.cell_relative_size_to_pixels(
    seq(1, 5, length.out = 100)
  )
  pixel_groups <- pixelatorR:::.cell_rgl_size_groups(pixel_sizes)
  expect_equal(anyNA(pixel_groups$keys), FALSE)
  expect_equal(max(pixel_groups$keys) <= 20L, TRUE)
  expect_equal(min(pixel_groups$keys[pixel_groups$keys > 0L]) >= 1L, TRUE)
  pixel_ids <- pixelatorR:::.cell_rgl_points(
    x = 1:100,
    y = 1:100,
    z = 1:100,
    color = rep("#000000", 100),
    size = pixel_sizes,
    alpha = rep(1, 100)
  )
  expect_equal(length(pixel_ids) <= 20L, TRUE)
  rgl::close3d()

  wide_categorical <- tibble::tibble(
    x = seq_len(20),
    y = 0,
    z = 0,
    cell = paste0("c", seq_len(20)),
    group = "a"
  ) |>
    cell_plot(color = group) |>
    cell_grid(cols = cell)
  legend_log <- tempfile()
  legend_con <- file(legend_log, open = "wt")
  sink(legend_con, type = "message")
  wide_device <- cell_plot_rgl(wide_categorical)
  sink(type = "message")
  close(legend_con)
  expect_equal(as.integer(wide_device), as.integer(rgl::cur3d()))
  expect_equal(
    any(grepl("figure margins too large", readLines(legend_log), fixed = TRUE)),
    FALSE
  )
  unlink(legend_log)
  rgl::close3d()

  flat_legend <- list(
    title = "marker",
    limits = c(2, 2),
    colors = c("#000000", "#FFFFFF")
  )
  rgl::open3d()
  expect_no_error({
    pixelatorR:::.cell_rgl_draw_colorbar(
      flat_legend,
      text_color = "black",
      text_size = 11
    )
    pixelatorR:::.cell_rgl_draw_colorbar(
      list(
        title = "marker",
        limits = c(NA_real_, NA_real_),
        colors = flat_legend$colors
      ),
      text_color = "black",
      text_size = 11
    )
  })
  expect_true(any(vapply(
    rgl::scene3d()$objects,
    function(object) identical(object$type, "quads"),
    logical(1)
  )))
  rgl::close3d()

  na_facet <- tibble::tibble(
    x = c(0, 1),
    y = c(0, 1),
    z = c(0, 1),
    panel = c("a", NA_character_),
    marker = c(1, 1)
  ) |>
    cell_plot(color = marker) |>
    cell_node_scale_color(colors = c("black", "white"), limits = c(1, 1)) |>
    cell_grid(cols = panel) |>
    cell_plot_rgl()
  expect_equal(na_facet, as.integer(rgl::cur3d()))
  expect_equal(rgl_plot_summary(na_facet)$n_panels, 2L)
})

test_that("rgl chrome is native text and legend geometry", {
  old_options <- options(rgl.useNULL = TRUE)
  on.exit(options(old_options), add = TRUE)
  on.exit(try(rgl::close3d(), silent = TRUE), add = TRUE)

  faceted <- tibble::tibble(
    x = c(0, 1, 2, 3),
    y = c(1, 0, 2, 3),
    z = c(-1, 1, 0, 2),
    marker = c(0, 2, 1, 3),
    cell = c("a", "b", "a", "b"),
    panel = c("m1", "m1", "m2", "m2")
  ) |>
    cell_plot(color = marker) |>
    cell_node_scale_color(colors = c("black", "white"), limits = c(0, 3)) |>
    cell_grid(cols = cell, rows = panel) |>
    cell_annotation(title = "Spectral layout") |>
    cell_plot_rgl()

  scene <- rgl::scene3d()
  subscene_texts <- function(subscene) {
    texts <- lapply(as.character(subscene$objects), function(id) {
      object <- scene$objects[[id]]
      if (is.null(object) || !identical(object$type, "text")) {
        return(character())
      }
      return(as.character(object$texts))
    })
    return(unlist(texts, use.names = FALSE))
  }
  subscenes <- scene$rootSubscene$subscenes
  texts <- unlist(lapply(subscenes, subscene_texts), use.names = FALSE)
  expect_true(all(c("a", "b", "m1", "m2", "Spectral layout") %in% texts))

  find_subscene <- function(label) {
    matches <- Filter(
      function(subscene) label %in% subscene_texts(subscene),
      subscenes
    )
    return(matches[[1]])
  }
  column_strip <- find_subscene("a")
  row_strip <- find_subscene("m1")
  expect_equal(
    column_strip$par3d$viewport[[4]],
    row_strip$par3d$viewport[[3]]
  )
  expect_equal(
    row_strip$par3d$userMatrix,
    rgl::rotationMatrix(
      pixelatorR:::.cell_plot_row_strip_angle * pi / 180,
      0,
      0,
      1
    )
  )
  expect_equal(column_strip$par3d$userMatrix, diag(4))

  legend <- Find(
    function(subscene) {
      types <- rgl_subscene_object_types(subscene, scene)
      return(any(types == "quads"))
    },
    subscenes
  )
  expect_false(is.null(legend))
  expect_equal(legend$par3d$mouseMode, c(
    none = "none",
    left = "none",
    right = "none",
    middle = "none",
    wheel = "none"
  ))
  expect_equal(as.integer(legend$par3d$listeners), as.integer(legend$id))
})

