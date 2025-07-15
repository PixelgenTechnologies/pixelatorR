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
    render_rotating_layout(xyz, gif_file, max_degree = 20)
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
})
