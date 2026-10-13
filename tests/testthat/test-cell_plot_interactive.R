plot_data <- tibble::tibble(
  x = c(0, 1),
  y = c(1, 0),
  z = c(-1, 1),
  marker = c(0, 2)
)
recipe <- cell_plot(plot_data, color = marker)

test_that("cell_plot_interactive dispatches to the requested renderer", {
  skip_if_not_installed("rgl")
  skip_if_not_installed("plotly")
  old_options <- options(rgl.useNULL = TRUE)
  on.exit(options(old_options), add = TRUE)

  default_widget <- cell_plot_interactive(recipe)
  rgl_widget <- cell_plot_interactive(recipe, renderer = "rgl")
  plotly_widget <- cell_plot_interactive(recipe, renderer = "plotly")

  expect_equal(
    list(
      default = class(default_widget),
      rgl = class(rgl_widget),
      plotly = class(plotly_widget)
    ),
    list(
      default = c("rglWebGL", "htmlwidget"),
      rgl = c("rglWebGL", "htmlwidget"),
      plotly = c("plotly", "htmlwidget")
    )
  )
})

test_that("cell_plot_interactive fails with invalid input", {
  expect_error(cell_plot_interactive(plot_data))
  expect_error(cell_plot_interactive(recipe, renderer = "ggplot"))
  expect_error(cell_plot_interactive(recipe, renderer = "plot"))
  expect_error(cell_plot_interactive(recipe, renderer = c("rgl", "ggplot")))
})
