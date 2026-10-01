plotly_trace_summary <- function(trace) {
  marker <- trace$marker
  return(list(
    type = trace$type,
    scene = trace$scene,
    x = as.numeric(unname(trace$x)),
    y = as.numeric(unname(trace$y)),
    z = as.numeric(unname(trace$z)),
    name = trace$name,
    legendgroup = trace$legendgroup,
    showlegend = isTRUE(trace$showlegend),
    visible = trace$visible,
    hoverinfo = trace$hoverinfo,
    marker = list(
      color = unname(as.character(marker$color)),
      size = as.numeric(unname(marker$size)),
      opacity = marker$opacity,
      showscale = isTRUE(marker$showscale),
      cmin = marker$cmin,
      cmax = marker$cmax,
      colorscale = marker$colorscale,
      colorbar_title = marker$colorbar$title$text
    )
  ))
}

plotly_scene_summary <- function(scene) {
  return(list(
    domain = scene$domain,
    range = list(
      x = scene$xaxis$range,
      y = scene$yaxis$range,
      z = scene$zaxis$range
    ),
    visible = list(
      x = scene$xaxis$visible,
      y = scene$yaxis$visible,
      z = scene$zaxis$visible
    ),
    aspectmode = scene$aspectmode,
    bgcolor = scene$bgcolor
  ))
}

plotly_plot_summary <- function(plot) {
  built <- plotly::plotly_build(plot)
  scene_names <- grep("^scene[0-9]*$", names(built$x$layout), value = TRUE)
  return(list(
    n_traces = length(built$x$data),
    traces = lapply(built$x$data, plotly_trace_summary),
    scenes = stats::setNames(
      lapply(scene_names, function(name) {
        plotly_scene_summary(built$x$layout[[name]])
      }),
      scene_names
    ),
    annotations = lapply(built$x$layout$annotations, function(a) a$text),
    title = built$x$layout$title$text,
    paper_bgcolor = built$x$layout$paper_bgcolor,
    font = built$x$layout$font[c("color", "size")],
    showlegend = built$x$layout$showlegend,
    legend_title = built$x$layout$legend$title$text
  ))
}

test_that("Interactive cell plots work as expected", {
  skip_if_not_installed("plotly")
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
    cell_plot_interactive()

  expect_s3_class(continuous, "plotly")
  expect_equal(
    plotly_plot_summary(continuous),
    list(
      n_traces = 2L,
      traces = list(
        list(
          type = "scatter3d",
          scene = "scene",
          x = c(0, 1),
          y = c(1, 0),
          z = c(-1, 1),
          name = NULL,
          legendgroup = NULL,
          showlegend = FALSE,
          visible = NULL,
          hoverinfo = NULL,
          marker = list(
            color = c("rgba(0,0,0,0.36)", "rgba(255,255,255,1)"),
            size = c(2, 6) * 96 / 25.4,
            opacity = NULL,
            showscale = FALSE,
            cmin = NULL,
            cmax = NULL,
            colorscale = NULL,
            colorbar_title = NULL
          )
        ),
        list(
          type = "scatter3d",
          scene = "scene",
          x = -0.05,
          y = -0.05,
          z = -1.1,
          name = NULL,
          legendgroup = NULL,
          showlegend = FALSE,
          visible = NULL,
          hoverinfo = "skip",
          marker = list(
            color = "0",
            size = 0,
            opacity = 0,
            showscale = TRUE,
            cmin = 0,
            cmax = 2,
            colorscale = list(list(0, "#000000"), list(1, "#FFFFFF")),
            colorbar_title = "marker"
          )
        )
      ),
      scenes = list(
        scene = list(
          domain = list(x = c(0, 1), y = c(0, 1)),
          range = list(
            x = c(-0.05, 1.05),
            y = c(-0.05, 1.05),
            z = c(-1.1, 1.1)
          ),
          visible = list(x = FALSE, y = FALSE, z = FALSE),
          aspectmode = "data",
          bgcolor = "navy"
        )
      ),
      annotations = list(),
      title = "Cells<br><sup>demo</sup>",
      paper_bgcolor = "navy",
      font = list(color = "white", size = 12),
      showlegend = FALSE,
      legend_title = "marker"
    )
  )

  constant <- cell_plot(plot_data) |>
    cell_plot_interactive()
  constant_marker <- plotly_plot_summary(constant)$traces[[1]]$marker
  expect_equal(
    constant_marker[c("color", "size")],
    list(
      color = c("rgba(229,229,229,1)", "rgba(229,229,229,1)"),
      size = c(1, 1) * 96 / 25.4
    )
  )

  mixed <- cell_plot(plot_data, color = cell) |>
    cell_node_scale_color(colors = c(a = "red", b = "blue")) |>
    cell_plot_interactive()
  expect_equal(
    plotly_plot_summary(mixed)$traces,
    list(
      list(
        type = "scatter3d",
        scene = "scene",
        x = c(0, 1),
        y = c(1, 0),
        z = c(-1, 1),
        name = NULL,
        legendgroup = NULL,
        showlegend = FALSE,
        visible = NULL,
        hoverinfo = NULL,
        marker = list(
          color = c("rgba(0,0,255,1)", "rgba(255,0,0,1)"),
          size = c(1, 1) * 96 / 25.4,
          opacity = NULL,
          showscale = FALSE,
          cmin = NULL,
          cmax = NULL,
          colorscale = NULL,
          colorbar_title = NULL
        )
      ),
      list(
        type = "scatter3d",
        scene = "scene",
        x = -0.05,
        y = -0.05,
        z = -1.1,
        name = "a",
        legendgroup = "a",
        showlegend = TRUE,
        visible = "legendonly",
        hoverinfo = "skip",
        marker = list(
          color = "#FF0000",
          size = 8,
          opacity = NULL,
          showscale = FALSE,
          cmin = NULL,
          cmax = NULL,
          colorscale = NULL,
          colorbar_title = NULL
        )
      ),
      list(
        type = "scatter3d",
        scene = "scene",
        x = -0.05,
        y = -0.05,
        z = -1.1,
        name = "b",
        legendgroup = "b",
        showlegend = TRUE,
        visible = "legendonly",
        hoverinfo = "skip",
        marker = list(
          color = "#0000FF",
          size = 8,
          opacity = NULL,
          showscale = FALSE,
          cmin = NULL,
          cmax = NULL,
          colorscale = NULL,
          colorbar_title = NULL
        )
      )
    )
  )

  categorical <- cell_plot(plot_data, color = cell) |>
    cell_node_scale_color(colors = c(a = "red", b = "blue")) |>
    cell_grid(cols = cell) |>
    cell_plot_interactive()
  expect_equal(
    plotly_plot_summary(categorical),
    list(
      n_traces = 4L,
      traces = list(
        list(
          type = "scatter3d",
          scene = "scene",
          x = 1,
          y = 0,
          z = 1,
          name = NULL,
          legendgroup = NULL,
          showlegend = FALSE,
          visible = NULL,
          hoverinfo = NULL,
          marker = list(
            color = "rgba(255,0,0,1)",
            size = 1 * 96 / 25.4,
            opacity = NULL,
            showscale = FALSE,
            cmin = NULL,
            cmax = NULL,
            colorscale = NULL,
            colorbar_title = NULL
          )
        ),
        list(
          type = "scatter3d",
          scene = "scene",
          x = -0.05,
          y = -0.05,
          z = -1.1,
          name = "a",
          legendgroup = "a",
          showlegend = TRUE,
          visible = "legendonly",
          hoverinfo = "skip",
          marker = list(
            color = "#FF0000",
            size = 8,
            opacity = NULL,
            showscale = FALSE,
            cmin = NULL,
            cmax = NULL,
            colorscale = NULL,
            colorbar_title = NULL
          )
        ),
        list(
          type = "scatter3d",
          scene = "scene",
          x = -0.05,
          y = -0.05,
          z = -1.1,
          name = "b",
          legendgroup = "b",
          showlegend = TRUE,
          visible = "legendonly",
          hoverinfo = "skip",
          marker = list(
            color = "#0000FF",
            size = 8,
            opacity = NULL,
            showscale = FALSE,
            cmin = NULL,
            cmax = NULL,
            colorscale = NULL,
            colorbar_title = NULL
          )
        ),
        list(
          type = "scatter3d",
          scene = "scene2",
          x = 0,
          y = 1,
          z = -1,
          name = NULL,
          legendgroup = NULL,
          showlegend = FALSE,
          visible = NULL,
          hoverinfo = NULL,
          marker = list(
            color = "rgba(0,0,255,1)",
            size = 1 * 96 / 25.4,
            opacity = NULL,
            showscale = FALSE,
            cmin = NULL,
            cmax = NULL,
            colorscale = NULL,
            colorbar_title = NULL
          )
        )
      ),
      scenes = list(
        scene = list(
          domain = list(x = c(0.02, 0.48), y = c(0.04, 0.96)),
          range = list(
            x = c(-0.05, 1.05),
            y = c(-0.05, 1.05),
            z = c(-1.1, 1.1)
          ),
          visible = list(x = FALSE, y = FALSE, z = FALSE),
          aspectmode = "data",
          bgcolor = "white"
        ),
        scene2 = list(
          domain = list(x = c(0.52, 0.98), y = c(0.04, 0.96)),
          range = list(
            x = c(-0.05, 1.05),
            y = c(-0.05, 1.05),
            z = c(-1.1, 1.1)
          ),
          visible = list(x = FALSE, y = FALSE, z = FALSE),
          aspectmode = "data",
          bgcolor = "white"
        )
      ),
      annotations = list("a", "b"),
      title = NULL,
      paper_bgcolor = "white",
      font = list(color = "black", size = 11),
      showlegend = TRUE,
      legend_title = "cell"
    )
  )

  expect_error(cell_plot(dplyr::select(plot_data, -z)))

  renamed_legend <- cell_plot(plot_data, color = cell) |>
    cell_annotation(legend_title = "Cell ID") |>
    cell_plot_interactive()
  renamed_summary <- plotly_plot_summary(renamed_legend)
  expect_equal(renamed_summary$legend_title, "Cell ID")
  expect_equal(
    renamed_summary$traces[[1]]$marker$colorbar_title,
    NULL
  )
})

test_that("Interactive cell plots ignore depth sizing and arrangement", {
  skip_if_not_installed("plotly")
  plot_data <- tibble::tibble(
    x = c(1, 0),
    y = c(1, 0),
    z = c(-1, 1)
  )

  without_depth <- cell_plot(plot_data, size = 2, depth = NULL) |>
    cell_plot_interactive()
  with_depth <- cell_plot(plot_data, size = 2) |>
    cell_node_depth(focal_distance = 0.2) |>
    cell_plot_interactive()
  arranged <- cell_plot(plot_data, arrange = x) |>
    cell_plot_interactive()

  expect_equal(
    plotly_plot_summary(without_depth)$traces[[1]]$marker$size,
    c(2, 2) * 96 / 25.4
  )
  expect_equal(
    plotly_plot_summary(with_depth)$traces[[1]]$marker$size,
    plotly_plot_summary(without_depth)$traces[[1]]$marker$size
  )
  expect_equal(
    plotly_plot_summary(arranged)$traces[[1]]$x,
    c(1, 0)
  )
})

test_that("Interactive cell plot grids work as expected", {
  skip_if_not_installed("plotly")
  na_factor <- tibble::tibble(
    x = 1:3,
    y = 1:3,
    z = 1:3,
    panel = factor(c("b", NA, "a"), levels = c("a", "b", "unused"))
  ) |>
    cell_plot() |>
    cell_grid(rows = panel) |>
    cell_plot_interactive()

  na_summary <- plotly_plot_summary(na_factor)
  expect_equal(
    na_summary$annotations,
    list("a", "b", "unused", "NA")
  )
  expect_equal(
    lapply(na_summary$traces, function(trace) {
      list(scene = trace$scene, x = trace$x, hoverinfo = trace$hoverinfo)
    }),
    list(
      list(scene = "scene", x = 3, hoverinfo = NULL),
      list(scene = "scene2", x = 1, hoverinfo = NULL),
      list(scene = "scene3", x = 2, hoverinfo = "skip"),
      list(scene = "scene4", x = 2, hoverinfo = NULL)
    )
  )
  expect_equal(
    unname(vapply(na_summary$scenes, function(scene) scene$aspectmode, character(1))),
    rep("data", 4)
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
      cell_plot_interactive()
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
      cell_plot_interactive()
  )
})
