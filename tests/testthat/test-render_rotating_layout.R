set.seed(123)
cg <- (ReadPNA_Seurat(minimal_pna_pxl_file(), verbose = FALSE) %>%
  LoadCellGraphs(cells = colnames(.)[1], add_layouts = TRUE, verbose = FALSE) %>%
  CellGraphs())[[1]]
xyz <- cg@layout$wpmds_3d
xyz$node_val <- runif(nrow(xyz))
gif_file <- fs::file_temp(ext = "gif")
mp4_file <- fs::file_temp(ext = "mp4")

test_that("render_rotating_layout works as expected", {
  # Simple case gif
  expect_no_error(render_rotating_layout(xyz, gif_file, frames = 2, show_first_frame = FALSE))

  # Simple case mp4
  expect_no_error(render_rotating_layout(xyz, mp4_file, frames = 2, show_first_frame = FALSE))

  # Change colors
  expect_no_error(render_rotating_layout(xyz, gif_file,
    frames = 2, show_first_frame = FALSE,
    colors = c("red", "blue")
  ))
  # Change opacity
  expect_no_error(render_rotating_layout(xyz, gif_file,
    frames = 2, show_first_frame = FALSE,
    pt_opacity = 0.1
  ))

  # Change point size
  expect_no_error(render_rotating_layout(xyz, gif_file,
    frames = 2, show_first_frame = FALSE,
    pt_size = 2
  ))

  # Change max_degree
  expect_no_error(render_rotating_layout(xyz, gif_file,
    frames = 2, show_first_frame = FALSE,
    max_degree = 120
  ))

  # Change center_zero
  expect_no_error(render_rotating_layout(xyz, gif_file,
    frames = 2, show_first_frame = FALSE,
    center_zero = TRUE
  ))

  # Change padding
  expect_no_error(render_rotating_layout(xyz, gif_file,
    frames = 2, show_first_frame = FALSE,
    pad = 0.5
  ))

  # Change width / height
  expect_no_error(render_rotating_layout(xyz, gif_file,
    frames = 2, show_first_frame = FALSE,
    width = 300, height = 300
  ))
  gif <- magick::image_read(gif_file)
  magick::image_info(gif)$width %>% expect_equal(300)
  magick::image_info(gif)$height %>% expect_equal(300)

  # Change resolution
  expect_no_error(render_rotating_layout(xyz, gif_file,
    frames = 2, show_first_frame = FALSE,
    res = 100, width = 300, height = 300
  ))

  # Use base R graphics
  expect_no_error(render_rotating_layout(xyz, gif_file,
    frames = 2, show_first_frame = FALSE,
    graphics_use = "base"
  ))

  # Use base R graphics
  expect_no_error(render_rotating_layout(xyz, gif_file,
    frames = 2, show_first_frame = FALSE,
    ggplot_theme = ggplot2::theme_classic()
  ))

  # Use illumination mask with ggplot2
  expect_no_error(render_rotating_layout(xyz, gif_file,
    frames = 2, show_first_frame = FALSE,
    use_illumination = TRUE
  ))

  # Use illumination without normalizing the mask
  expect_no_error(render_rotating_layout(xyz, gif_file,
    frames = 2, show_first_frame = FALSE,
    use_illumination = TRUE,
    normalize_illumination = FALSE
  ))

  # Use illumination mask with base R graphics
  expect_no_error(render_rotating_layout(xyz, gif_file,
    frames = 2, show_first_frame = FALSE,
    graphics_use = "base",
    use_illumination = TRUE
  ))

  # Use illumination mask with factor node_val and shadow palette
  xyz_factor <- xyz %>%
    mutate(node_val = factor(ifelse(node_val > 0.9, "high", "low")))
  expect_no_error(render_rotating_layout(xyz_factor, gif_file,
    frames = 2, show_first_frame = FALSE,
    graphics_use = "ggplot2",
    colors = c("red", "gray95"),
    use_illumination = TRUE,
    illumination_shadow_colors = PixelgenGradient(100, "NaturalBlue")
  ))

  # Use shadow palette blending with ggplot2
  expect_no_error(render_rotating_layout(xyz, gif_file,
    frames = 2, show_first_frame = FALSE,
    use_illumination = TRUE,
    illumination_shadow_colors = c("#1F385A", "#C1DBE0")
  ))

  # Use shadow palette blending with base R graphics
  expect_no_error(render_rotating_layout(xyz, gif_file,
    frames = 2, show_first_frame = FALSE,
    graphics_use = "base",
    use_illumination = TRUE,
    illumination_shadow_colors = c("#1F385A", "#C1DBE0")
  ))
})


test_that("render_rotating_layout fails with invalid input", {
  # Invalid cell_col
  expect_error(
    render_rotating_layout(xyz, gif_file, cell_col = "Invalid")
  )

  # Invalid marker_col
  expect_error(
    render_rotating_layout(xyz, gif_file, marker_col = "Invalid")
  )

  # Invalid pt_opacity
  expect_error(
    render_rotating_layout(xyz, gif_file, pt_opacity = -1)
  )

  # Invalid pt_size
  expect_error(
    render_rotating_layout(xyz, gif_file, pt_size = -1)
  )

  # Invalid colors
  expect_error(
    render_rotating_layout(xyz, gif_file, colors = c("Invalid"))
  )

  # Invalid max_degree
  expect_error(
    render_rotating_layout(xyz, gif_file, max_degree = -20)
  )

  # Invalid center_zero
  expect_error(
    render_rotating_layout(xyz, gif_file, center_zero = "Invalid")
  )

  # Invalid frames
  expect_error(
    render_rotating_layout(xyz, gif_file, frames = -1)
  )

  # Invalid pad
  expect_error(
    render_rotating_layout(xyz, gif_file, pad = -1)
  )

  # Invalid show_first_frame
  expect_error(
    render_rotating_layout(xyz, gif_file, show_first_frame = "Invalid")
  )

  # Invalid width
  expect_error(
    render_rotating_layout(xyz, gif_file, width = -1)
  )

  # Invalid height
  expect_error(
    render_rotating_layout(xyz, gif_file, height = -1)
  )

  # Invalid res
  expect_error(
    render_rotating_layout(xyz, gif_file, res = "Invalid")
  )

  # Invalid delay
  expect_error(
    render_rotating_layout(xyz, gif_file, delay = -1)
  )

  # Invalid ggplot_theme
  expect_error(
    render_rotating_layout(xyz, gif_file, ggplot_theme = "Invalid")
  )

  # Invalid title
  expect_error(
    render_rotating_layout(xyz, gif_file, title = NULL)
  )

  # Invalid bg
  expect_error(
    render_rotating_layout(xyz, gif_file, bg = "Invalid")
  )

  # Invalid label_grid_axes
  expect_error(render_rotating_layout(xyz, gif_file, label_grid_axes = "Invalid"))

  # Invalid graphics_use
  expect_error(
    render_rotating_layout(xyz, gif_file, graphics_use = "Invalid")
  )

  # Invalid illumination_ambient
  expect_error(
    render_rotating_layout(xyz, gif_file, illumination_ambient = -1)
  )

  # Invalid illumination_sat_boost
  expect_error(
    render_rotating_layout(xyz, gif_file, illumination_sat_boost = -1)
  )

  # Invalid illumination_shadow_colors
  expect_error(
    render_rotating_layout(
      xyz,
      gif_file,
      use_illumination = TRUE,
      illumination_shadow_colors = c("Invalid")
    )
  )
})

test_that(".apply_palette_illumination blends colors as expected", {
  base <- "#FF0000"
  shadow <- "#0000FF"

  full_light <- pixelatorR:::.apply_palette_illumination(
    base, shadow,
    illum_mask = 1, ambient_intensity = 0
  )
  expect_equal(full_light, base)

  full_shadow <- pixelatorR:::.apply_palette_illumination(
    base, shadow,
    illum_mask = 0, ambient_intensity = 0
  )
  expect_equal(full_shadow, shadow)

  expect_equal(
    pixelatorR:::.apply_palette_illumination(
      rep(base, 9), shadow,
      illum_mask = seq(0, 1, length.out = 9), ambient_intensity = 0
    ),
    c(
      "#0000FF", "#1F00DF", "#3F00BF", "#5F009F", "#7F007F", "#9F005F",
      "#BF003F", "#DF001F", "#FF0000"
    )
  )
})

test_that("illumination helpers validate inputs", {
  expect_error(
    pixelatorR:::.apply_hsv_illumination("Invalid", 0.5),
    "valid color|color"
  )
  expect_error(
    pixelatorR:::.apply_hsv_illumination("#FF0000", c(0.5, 0.6)),
    "same length|length"
  )
  expect_error(
    pixelatorR:::.apply_palette_illumination("#FF0000", "#0000FF", 1.5),
    "between|range"
  )
  expect_error(
    pixelatorR:::.apply_palette_illumination(
      c("#FF0000", "#00FF00"),
      c("#0000FF", "#0000FF", "#0000FF"),
      c(0.5, 0.5)
    ),
    "same length|length"
  )
})

test_that(".node_val_lims works as expected", {
  expect_equal(
    pixelatorR:::.node_val_lims(c(-2, 5), center_zero = FALSE),
    c(-2, 5)
  )
  expect_equal(
    pixelatorR:::.node_val_lims(
      c(-2, 5),
      marker_limits = list(marker = c(-3, 4)),
      marker_id = "marker",
      center_zero = FALSE
    ),
    c(-3, 4)
  )
  expect_equal(
    pixelatorR:::.node_val_lims(c(-2, 5), center_zero = TRUE),
    c(-5, 5)
  )

  expect_equal(nrow(pixelatorR:::.legend_carrier_data(runif(1e4))), 2)
  expect_equal(
    nrow(pixelatorR:::.legend_carrier_data(factor(rep(c("a", "b"), 5e3)))),
    2
  )
})

test_that(".node_val_color_scale uses visible legend keys with illumination", {
  scale <- pixelatorR:::.node_val_color_scale(
    factor(c("high", "low")),
    colors = c("red", "blue"),
    center_zero = FALSE,
    use_illumination = TRUE,
    pt_size = 2
  )
  guide <- scale$guide
  expect_s3_class(guide, "GuideLegend")
  expect_equal(scale$guide$params$override.aes$alpha, 1)
  expect_equal(scale$guide$params$override.aes$size, 2)
})

test_that(".compute_illuminated_point_colors works as expected", {
  # Numeric node_val
  df <- tibble(
    node_val = c(-1, 0, 1),
    illumination = c(0, 0.5, 1),
    marker = "marker"
  )

  marker_limits <-
    list(marker = tibble(min = -1, max = 1))

  expect_no_error(
    adjusted_colors <-
      .compute_illuminated_point_colors(
        df,
        colors = c("red", "blue"),
        marker_limits = marker_limits,
        marker_id = "marker",
        center_zero = TRUE,
        illumination_ambient = 0.5,
        illumination_sat_boost = 0.5
      )
  )

  expect_equal(
    adjusted_colors,
    c("#800000", "#980066", "#0000FF")
  )

  # Factor node_val
  df <- tibble(
    node_val = factor(c("val1", "val2", "val2", "val3")),
    illumination = c(0, 0.5, 0.75, 1),
    marker = "marker"
  )

  marker_limits <-
    list(marker = tibble(min = -1, max = 1))

  expect_no_error(
    adjusted_colors <-
      .compute_illuminated_point_colors(
        df,
        colors = c(val1 = "red", val3 = "blue", val2 = "white"),
        marker_limits = marker_limits,
        marker_id = "marker",
        center_zero = TRUE,
        illumination_ambient = 0.5,
        illumination_sat_boost = 0.5
      )
  )

  expect_equal(
    adjusted_colors,
    c("#800000", "#BFBFBF", "#DFDFDF", "#0000FF")
  )
})
