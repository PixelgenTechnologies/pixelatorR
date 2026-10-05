#' Render a cell plot as an interactive 3D widget
#'
#' Builds a [cell_plot()] recipe and draws it as an interactive 3D scatter
#' with one of two rendering backends, returned as an htmlwidget. Drag to
#' rotate and scroll to zoom. Markers use the rendered colors resolved by
#' [build_cell_plot()], including illumination when requested, while legends
#' are built from the unilluminated color scale metadata. Node sizes are
#' converted from relative units to marker diameters in pixels.
#'
#' Occlusion follows the scene camera, so markers closer to the current
#' viewpoint appear in front. The `arrange` and `depth` mappings and
#' [cell_coord_rotate()] are ignored. Panel grids support at most 10 rows and
#' 20 columns.
#'
#' `renderer = "rgl"` draws the scene with WebGL through [rgl]. The widget
#' has no fixed size: it fills the IDE viewer, a Quarto or R Markdown page,
#' or any other container and follows that container when it is resized.
#' Panels in a grid share one camera and move together. Titles, facet strips,
#' and legends are HTML drawn over the WebGL canvas, so they stay sharp at any
#' resolution and their text can be selected. Their size follows the theme
#' text size in points. Numeric color mappings get a colorbar and categorical
#' mappings a discrete legend. Continuous sizes are grouped into at most 20
#' size bins. Hover labels and interactive legend filtering are not available.
#'
#' `renderer = "plotly"` draws the scene with [plotly]. Each panel is a
#' separate scene with its own camera, and hovering a marker shows its
#' coordinates. Titles, facet strips, and legends are part of the Plotly
#' layout.
#'
#' @param object A `cell_plot` recipe.
#' @param renderer Rendering backend, `"rgl"` (default) or `"plotly"`. The
#' matching package must be installed.
#'
#' @return An htmlwidget. With `renderer = "rgl"` an rgl widget, with
#' `renderer = "plotly"` a Plotly widget.
#'
#' @seealso [cell_plot()], [print.cell_plot()], [cell_plot_animate()]
#'
#' @examplesIf interactive()
#' se <- ReadPNA_Seurat(minimal_pna_pxl_file())
#' se <- LoadCellGraphs(se, cells = colnames(se)[4], verbose = FALSE) |>
#'   ComputeLayout(layout_method = "spectral")
#'
#' cell_graph <- CellGraphs(se)[[4]]
#'
#' layout_data <- FetchLayoutData(cell_graph, vars = "CD82", layout_method = "spectral_3d")
#'
#' # WebGL scene with rgl
#' cell_plot(layout_data, color = CD82) |>
#'   cell_plot_interactive()
#'
#' # The same recipe drawn with Plotly
#' cell_plot(layout_data, color = CD82) |>
#'   cell_plot_interactive(renderer = "plotly")
#'
#' @export
cell_plot_interactive <- function(object, renderer = c("rgl", "plotly")) {
  .validate_cell_plot(object)
  renderer <- rlang::arg_match(renderer)
  switch(renderer,
    rgl = expect_rgl(),
    plotly = expect_plotly()
  )

  object$mapping$arrange <- NULL
  built <- build_cell_plot(object)
  widget <- switch(renderer,
    rgl = .cell_rgl_html_widget(built),
    plotly = .render_cell_plot_plotly(built)
  )
  return(widget)
}
