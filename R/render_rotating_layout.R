#' Create a rotating 3D layout video
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' \code{render_rotating_layout} can be used to generate a rotating 3D scatter
#' plot from a tibble with node layout coordinates. \code{render_rotating_layout}
#' offers a number of customization options, including the option to facet the
#' plot by additional columns in the tibble. See the sections below for more details
#' on how to adjust the appearance.
#'
#' @section cells and markers:
#' The \code{cell_col} and \code{marker_col} arguments can be used to facet
#' the plot by additional columns in the tibble. For instance, if the tibble
#' contains layout coordinates for multiple cells and markers, two additional
#' columns should specify these groupings. In the final plot, each cell will
#' plotted in its own coordinate system to make sure that the layout fills
#' the plot area. Similarly, each marker will share the same color scale to make
#' sure that the colors are comparable across cells.
#'
#' @section graphic parameters:
#' The appearance of the plot can be customized using a number of graphics
#' parameters. Increasing the number of \code{frames} while decreasing the
#' delay will create a smoother transition but will take much longer to render.
#' The \code{width} and \code{height} parameters should be increased manually
#' to accommodate more cells and markers or to increase resolution.
#' The \code{width} and \code{height} parameters should be adjusted to match
#' the selected dpi (\code{res}). Note that when these parameters are increased,
#' the rendering time can increase drastically. Unfortunately, drawing ggplot
#' objects to the graphics device is relatively slow so the rendering might take
#' several minutes. Using base R graphics instead of ggplot2 can speed up the
#' rendering significantly, but it will limit the customization options.
#'
#' Since the rendering is expensive, we might want to inspect the first frame to
#' make sure that the plot looks as expected. If the \code{show_first_frame} argument
#' is set to \code{TRUE}, which is the default, the first frame will be displayed
#' in the viewer before proceeding with the rendering. In interactive sessions,
#' you will be prompted to enter 'Yes' to continue if you are satisfied with the
#' first frame.
#'
#' @section parallel processing:
#' To speed up the rendering, the \code{cl} argument can be used to specify a
#' number of threads. This feature only works on unix-based systems.
#' This can significantly reduce the rendering time; however, it might create
#' artefacts in the output video. The default is \code{NULL}, which corresponds
#' to a single-threaded operation (sequential processing).
#'
#' @section graphics_use:
#' The \code{graphics_use} argument can be used to switch between the default
#' "base" and the "ggplot2" graphics system. The "base" option is much faster
#' but has fewer customization options. For instance, it is currently not possible to
#' get row and column labels with the "base" graphics option.
#'
#' @param data A tibble (\code{tbl_df}) with columns 'x', 'y', 'z',
#' and 'node_val'. The 'node_val' column can be either a numeric or a
#' factor.
#' @param file A character string specifying the path to the output file.
#' The video format will be determined by the file extension. If the extension is
#' '.gif', the \code{gifski} R package will be used to render a gif. Other file
#' formats (such as mp4, mkv, mov, or flv) will be rendered using the \code{av}
#' R package. Make sure to use a file format supported by either \code{gifski}
#' or \code{av}.
#' @param cell_col,marker_col A character string specifying columns to facet the
#' plot by. \code{cell_col} should be a column with unique cell identifiers
#' and \code{marker_col} should be a column with unique protein identifiers.
#' @param pt_opacity A numeric value between 0 and 1 specifying the opacity
#' of the points.
#' @param pt_size A numeric value indicating the maximum size of the points.
#' @param colors A vector of valid colors. If 'node_val' is a factor,
#' the length of 'colors' should be equal to the number of unique levels.
#' If 'node_val' is numeric, the colors will be used to create a gradient
#' color scale.
#' @param max_degree A numeric value between 90 and 360. The maximum angle
#' of rotation around the z axis. Default is 360, which corresponds to a
#' full turn around the z axis.
#' @param center_zero A logical value indicating whether the color gradient
#' should be centered around zero.
#' @param scale_layout A logical value indicating whether the layouts should be
#' such that they all fit the same bounding box. This is typically necessary
#' as layouts vary significantly in size.
#' @param frames A positive numeric value indicating the number of frames
#' to render. More frames will result in a smoother animation but will
#' increase the size of the output file and rendering time.
#' @param pad A numeric value between 0 and 1. The amount of padding to add
#' to the axis range. Default is 0.1 which corresponds to an expansion of 10%
#' along all three axes.
#' @param width,height An integer specifying the width and height of the
#' video in pixels. Default is 500x500 which is suitable for a low resolution
#' GIF file.
#' @param delay A numeric value indicating the delay between frames in seconds.
#' Default is 1/20, which corresponds to 20 frames per second. The value should
#' be between 1/100 and 1.
#' @param show_first_frame A logical value indicating whether the first frame
#' should be displayed in the viewer before proceeding with the rendering. In
#' interactive sessions, you will be prompted to enter 'Yes' to continue if you
#' are satisfied with the first frame. Default is TRUE.
#' @param res A numeric value indicating the resolution of the png files in dpi.
#' See \code{\link[grDevices]{png}} for more details.
#' @param ggplot_theme A ggplot2 theme object. Default is \code{NULL}, which
#' corresponds to the default plot theme.
#' @param title A character string specifying the title of the plot.
#' @param bg A valid color string specifying the background color of the plot.
#' @param label_grid_axes A logical value indicating whether the rows and columns of
#' the plot grid should be labeled.
#' @param margin_widths A numeric vector of length 2 specifying the width of the
#' left and top margins relative to the width and height of the plot. Default is
#' c(0.1, 0.1), meaning that the left and top margins will take up 10% of the
#' width and height respectively.
#' @param use_facet_grid If set to TRUE, the plot will be faceted using
#' \code{\link[ggplot2]{facet_grid}}. Note that with this option, the node
#' colors will be mapped to a single color scale. If these selected markers
#' have very different dynamic ranges, it might be difficult to observe trends
#' in lowly expressed markers.
#' @param flip A logical value indicating whether the plot should be flipped
#' such that markers are arranged in columns and cells in rows.
#' @param graphics_use One of "ggplot2" or "base". The default, "base", is much
#' faster but has fewer customization options. The "ggplot2" option
#' comes with a number of customization options.
#' @param boomerang A logical value specifying whether the rotation should
#' be reversed after reaching the maximum degree to create a "boomerang" effect.
#' The frames are simply reversed and appended to the original frames.
#' The resulting video will have twice as many frames as specified by the
#' \code{frames} argument and keep the same frame rate. Note that if the output
#' is a gif file, it will double in size.
#' @param cl A number of threads or a cluster object created by
#' \code{\link[parallel]{makeCluster}}. The default is \code{NULL}, which
#' corresponds to a single-threaded operation (sequential processing).
#' @param keep_frames A logical value indicating whether the png files should
#' be kept after rendering the video. This option is useful if you want to
#' access to the individual frames later for further processing.
#'
#' @returns Exports an animation of a rotating 3D scatter plot.
#'
#' @examples
#' library(dplyr)
#' pxl_file <- minimal_pna_pxl_file()
#' se <- ReadPNA_Seurat(pxl_file)
#' se <- se %>%
#'   LoadCellGraphs(add_layouts = TRUE)
#'
#' # Create a gif from a 3D layout
#' cg <- CellGraphs(se)[[3]]
#' df <- cg@layout$wpmds_3d %>%
#'   mutate(node_val = cg@counts[, "CD3e"])
#' temp_gif <- fs::file_temp(ext = ".gif")
#' render_rotating_layout(
#'   data = df,
#'   file = temp_gif,
#'   colors = c("lightgrey", "red"),
#'   pt_size = 0.5,
#'   max_degree = 90,
#'   frames = 20,
#'   delay = 1 / 5,
#'   bg = "transparent",
#'   label_grid_axes = FALSE,
#'   show_first_frame = FALSE
#' )
#' magick::image_read(temp_gif)
#'
#' \dontrun{
#' # Include multiple facets
#' markers <- paste0("Marker", 1:3)
#' df <- lapply(1:2, function(i) {
#'   set.seed(i)
#'   cg <- simulate_bipartite_graph(n_nodes = 5e3, epsilon = 8) %>%
#'     add_binary_marker_counts_pol(marker_polarize = "Marker1", epsilon = 12)
#'   xyz <- layout_with_weighted_pmds(cg@cellgraph, dim = 3) %>%
#'     as_tibble(.name_repair = ~ c("x", "y", "z"))
#'   lapply(colnames(cg@counts)[1:2], function(m) {
#'     xyz %>%
#'       mutate(node_val = cg@counts[, m] %>% as.character()) %>%
#'       mutate(marker = m)
#'   }) %>%
#'     bind_rows() %>%
#'     mutate(cell = paste0("Cell", i))
#' }) %>% bind_rows()
#'
#' temp_gif <- fs::file_temp(ext = ".gif")
#' render_rotating_layout(
#'   data = df,
#'   pt_size = 0.5,
#'   width = 740,
#'   height = 650,
#'   cell_col = "cell",
#'   marker_col = "marker",
#'   colors = viridis::viridis(n = 2),
#'   file = temp_gif,
#'   max_degree = 90,
#'   frames = 100,
#'   delay = 1 / 20,
#'   res = 150,
#'   bg = "grey",
#'   show_first_frame = FALSE,
#'   cl = 9
#' )
#' magick::image_read(temp_gif)
#'
#' # Use base R graphics
#' temp_gif <- fs::file_temp(ext = ".gif")
#' render_rotating_layout(
#'   data = df,
#'   pt_size = 0.5,
#'   pt_opacity = 1,
#'   width = 740,
#'   height = 650,
#'   cell_col = "cell",
#'   marker_col = "marker",
#'   colors = viridis::viridis(n = 2),
#'   file = temp_gif,
#'   max_degree = 360,
#'   frames = 500,
#'   delay = 1 / 20,
#'   res = 150,
#'   bg = "grey",
#'   show_first_frame = FALSE,
#'   graphics_use = "base",
#'   cl = 9
#' )
#' magick::image_read(temp_gif)
#' }
#'
#' @export
#'
render_rotating_layout <- function(
  data,
  file,
  cell_col = NULL,
  marker_col = NULL,
  pt_opacity = 0.7,
  pt_size = 1,
  colors = RColorBrewer::brewer.pal(9, "Blues"),
  max_degree = 360,
  center_zero = FALSE,
  scale_layout = TRUE,
  frames = 500,
  pad = 0.1,
  show_first_frame = TRUE,
  width = 500,
  height = 500,
  res = 150,
  delay = 1 / 20,
  ggplot_theme = NULL,
  title = "",
  bg = "white",
  label_grid_axes = TRUE,
  margin_widths = c(0.1, 0.1),
  use_facet_grid = FALSE,
  flip = FALSE,
  graphics_use = c("base", "ggplot2"),
  boomerang = FALSE,
  cl = NULL,
  keep_frames = FALSE
) {
  if (fs::path_ext(file) == "gif") {
    rlang::check_installed("gifski")
  } else {
    rlang::check_installed("av")
  }

  # ragg is slightly faster than png, but needs to be installed
  if (!requireNamespace("ragg")) {
    warn(glue("For faster rendering, install the 'ragg' package."))
    dev_png <- png
  } else {
    dev_png <- ragg::agg_png
  }

  graphics_use <- match.arg(graphics_use, choices = c("base", "ggplot2"))
  if (graphics_use == "base") {
    cli_alert_warning(
      glue(
        "The following parameters are not supported with graphics_use = 'base':\n",
        "ggplot_theme, label_grid_axes, title\n",
        "Titles, text annotations and color bar will be missing from the plot."
      )
    )
  }

  .validate_render_rotating_layout_input_params(
    data, marker_col, cell_col, pt_opacity,
    pt_size, colors, max_degree, center_zero, scale_layout,
    frames, pad, show_first_frame, width, height,
    res, delay, ggplot_theme, title, bg,
    label_grid_axes, margin_widths, use_facet_grid,
    flip, boomerang
  )

  # Set variables if NULL
  cell_col <- cell_col %||% "cell"
  marker_col <- marker_col %||% "marker"
  # Add missing columns to data required for downstream processing
  if (!cell_col %in% names(data)) {
    data <- data %>%
      mutate(!!sym(cell_col) := cell_col)
  }
  if (!marker_col %in% names(data)) {
    data <- data %>%
      mutate(!!sym(marker_col) := marker_col)
  }

  # Create a vector of angles for each frame
  angles <- seq(0, max_degree, length.out = frames)

  # Create chunks to process
  indices_split <- seq_along(angles) %>% split(ceiling(seq_along(angles) / 3))

  # Create temporary directory where the frames will be saved
  tmp_dir <- fs::file_temp() %>% stringr::str_replace("file", "dir")
  dir.create(tmp_dir)

  # Get cell axes ranges
  # First, we group the tibble by the cell column, then we
  # split the tibble into a list of tibbles, one for each cell.
  # We then rescale the x, y, and z coordinates of each cell
  # such that the maximum absolute value of the coordinates is 1.
  # This is to ensure that all cells will fit to the same plot area.
  xyz_grouped <- data %>%
    select(x, y, z, node_val, all_of(c(marker_col, cell_col))) %>%
    group_by(across(all_of(cell_col)))
  group_d <- xyz_grouped %>% group_keys()
  xyz_list <- xyz_grouped %>%
    group_split() %>%
    set_names(nm = if (nrow(group_d) == 1) "cell" else group_d[, 1, drop = TRUE])
  if (scale_layout) {
    xyz_list <- lapply(xyz_list, function(xyz) {
      xyz_scaled <- xyz %>%
        select(x, y, z) %>%
        scale_layout(percentile = 1)
      xyz[, c("x", "y", "z")] <- xyz_scaled
      return(xyz)
    })
  }

  # Set the axis limits
  padded_max_lim <- 1 + 2 * pad
  xyz_limits <- c(-padded_max_lim, padded_max_lim)

  if (inherits(data$node_val, "numeric")) {
    # Get marker value limits
    # First, we group the tibble by the marker column, then we
    # split the tibble into a list of tibbles, one for each marker.
    # We then summarize the minimum and maximum values of the node_val
    # column for each marker. This is to make sure that each marker
    # has the same fixed color scale across different cells.
    marker_limits <- data %>%
      group_by(!!sym(marker_col)) %>%
      summarize(
        min = min(node_val, na.rm = TRUE),
        max = max(node_val, na.rm = TRUE), .groups = "keep"
      )
    marker_limits <- marker_limits %>%
      group_split() %>%
      set_names(nm = marker_limits[, 1, drop = TRUE]) %>%
      lapply(function(x) {
        x %>%
          select(min, max) %>%
          as.numeric()
      })
  } else {
    marker_limits <- NULL
  }

  # Create a nested list of tibbles, one for each cell, and one for each marker
  xyz_list_nested <- lapply(xyz_list, function(xyz) {
    xyz_list_cell <- xyz %>% group_by(across(all_of(marker_col)))
    group_d <- xyz_list_cell %>% group_keys()
    xyz_list_cell <- xyz_list_cell %>%
      group_split() %>%
      set_names(nm = if (nrow(group_d) == 1) "marker" else group_d[, 1, drop = TRUE])
    return(xyz_list_cell)
  })

  # Show the first frame if requested
  if (show_first_frame && interactive()) {
    .render_first_frame(
      xyz_list_nested, marker_limits, pt_opacity, pt_size,
      marker_col, cell_col, colors, center_zero, xyz_limits,
      tmp_dir, width, height, res, ggplot_theme, title, bg,
      use_facet_grid, flip, label_grid_axes, margin_widths,
      graphics_use, dev_png
    )
  }

  cli::cli_h3("Rendering frames...")
  # Loop over the indices chunks and create the frames
  dump <- pbapply::pblapply(indices_split, function(indices) {
    for (i in indices) {
      # Apply rotation and arrange values by the y axis
      # such that the nodes are stacked from back to front
      xyz_list_cur <- lapply(xyz_list_nested, function(xyz_cell) {
        xyz_cell <- lapply(xyz_cell, function(xyz) {
          xyz_rotated <- .rotate_around_z(xyz %>% select(x, y, z), angles[i])
          xyz[, c("x", "y", "z")] <- xyz_rotated
          xyz <- xyz %>%
            arrange(pick(all_of("y")))
          return(xyz)
        })
        return(xyz_cell)
      })

      # Current frame PNG file
      file <- file.path(
        tmp_dir,
        paste0("plot_", stringr::str_pad(i, side = "left", width = 4, "0"), ".png")
      )

      # Create and export plot using ggplot2
      if (graphics_use == "ggplot2") {
        p <- .scatter_3d_plot(
          xyz_list_cur,
          cell_col,
          marker_col,
          marker_limits,
          colors,
          center_zero,
          xyz_limits,
          pt_opacity,
          pt_size,
          use_facet_grid,
          flip,
          label_grid_axes,
          margin_widths
        )
        # Add ggplot theme
        if (!is.null(ggplot_theme)) {
          p <- p & ggplot_theme
        }
        if (title != "") {
          p <- p + plot_annotation(title = title)
        }
        ggsave(file, p, width = width, height = height, dpi = res, device = dev_png, units = "px", bg = bg)
      }

      # Create and export plot using base R graphics
      if (graphics_use == "base") {
        dev_png(file, width = width, height = height, res = res, units = "px", bg = bg)
        .scatter_3d_plot_base_R(
          xyz_list_cur,
          marker_limits,
          colors,
          center_zero,
          xyz_limits,
          pt_opacity,
          pt_size,
          bg,
          flip
        )
        dev.off()
      }
    }
  }, cl = cl)

  # Create video
  cli::cli_h3("Creating video...")
  png_files <- fs::dir_ls(tmp_dir, glob = "*.png")
  if (boomerang) {
    png_files <- c(png_files, rev(png_files))
  }

  if (fs::path_ext(file) == "gif") {
    # Use gifski if the file extension is .gif
    out_path <- gifski::gifski(png_files, file, width = width, height = height, delay = delay)
  } else {
    # Use av if the file extension is not .gif
    try_message <- try(
      {
        av::av_encode_video(png_files, file, framerate = 1 / delay)
      },
      silent = TRUE
    )
    # If the video creation fails, the most likely reason is that
    # the format is not supported. In this case, we provide a warning
    # and return the path to the directory containing the frames.
    if (inherits(try_message, "try-error")) {
      cli::cli_alert_danger(
        glue(
          "{try_message}\n",
          "Failed to generate the video. ",
          "Make sure to use an appropriate file extension. \n",
          "The frames can be found in:\n {tmp_dir}"
        )
      )
      return(invisible(NULL))
    }
  }

  # Remove temporary directory to clean up
  if (!keep_frames) {
    try_message <- try(
      {
        fs::dir_delete(tmp_dir)
      },
      silent = TRUE
    )
    if (inherits(try_message, "try-error")) {
      warn(glue("Failed to remove temporary directory:\n {tmp_dir}"))
    }
  } else {
    cli::cli_alert_info(glue("The frames can be found in:\n {tmp_dir}"))
  }

  cli::cli_alert_success(glue("The video has been saved to:\n {file}"))
}


#' Draws a scatter plot with simulated depth
#'
#' @param xyz_list_nested A list of tibbles (\code{tbl_df}) with columns 'x',
#' y', 'z', and 'node_val'.
#' @param cell_col,marker_col A character string indicating the column name
#' of the cell IDs or marker names.
#' @param marker_limits A list with numeric vectors of length 2. The range
#' of the marker values.
#' @param colors A vector of valid colors.
#' @param center_zero A logical value indicating whether the plot should be
#' centered around zero.
#' @param xyz_limits A numeric vector of length 2. The range of the x, y, and z
#' axes.
#' @param pt_opacity A numeric value between 0 and 1 indicating the opacity of
#' the points.
#' @param pt_size A numeric value indicating the maximum size of the points.
#' @param use_facet_grid A logical value indicating whether to use \code{facet_grid}
#' to construct the layout.
#' @param label_grid_axes A logical value indicating whether the axes should
#' be labeled.
#' @param margin_widths A numeric vector of length 2 indicating the width of the
#' margins relative to the plot size.
#'
#' @noRd
#'
.scatter_3d_plot <- function(
  xyz_list_nested,
  cell_col,
  marker_col,
  marker_limits,
  colors,
  center_zero,
  xyz_limits,
  pt_opacity,
  pt_size,
  use_facet_grid = FALSE,
  flip = FALSE,
  label_grid_axes = TRUE,
  margin_widths = c(0.1, 0.1)
) {
  if (use_facet_grid) {
    xyz_collapsed <- lapply(xyz_list_nested, function(xyz_cell) {
      xyz_cell %>%
        unname() %>%
        bind_rows()
    }) %>%
      unname() %>%
      bind_rows()
    p <- ggplot(xyz_collapsed, aes(x, z, size = y, color = node_val)) +
      geom_point(alpha = pt_opacity) +
      scale_size(range = c(0, pt_size), limits = xyz_limits) +
      coord_fixed() +
      theme_void() +
      guides(size = "none") +
      scale_x_continuous(limits = xyz_limits, expand = expansion()) +
      scale_y_continuous(limits = xyz_limits, expand = expansion()) +
      {
        if (inherits(xyz_collapsed$node_val, "numeric")) {
          max_abs_val <- max(abs(xyz_collapsed$node_val))
          lims <- if (center_zero) c(-max_abs_val, max_abs_val) else c(0, max_abs_val)
          scale_color_gradientn(colours = colors, limits = lims)
        } else {
          scale_color_manual(values = colors)
        }
      } +
      {
        if (flip) {
          facet_grid(reformulate(marker_col, cell_col), switch = "y")
        } else {
          facet_grid(reformulate(cell_col, marker_col), switch = "y")
        }
      }
  }

  if (!use_facet_grid) {
    # Create a list of ggplot objects per cell
    cell_plots <- lapply(names(xyz_list_nested), function(cell_id) {
      xyz_list_cell <- xyz_list_nested[[cell_id]]
      # Create a list of ggplot objects per marker
      marker_plots <- lapply(names(xyz_list_cell), function(marker_id) {
        p <- ggplot(
          xyz_list_cell[[marker_id]] %>%
            mutate(row_text = marker_id),
          aes(x, z, size = y, color = node_val)
        ) +
          geom_point(alpha = pt_opacity) +
          scale_size(range = c(0, pt_size), limits = xyz_limits) +
          coord_fixed() +
          theme_void() +
          guides(size = "none") +
          scale_x_continuous(limits = xyz_limits, expand = expansion()) +
          scale_y_continuous(limits = xyz_limits, expand = expansion()) +
          {
            if (inherits(xyz_list_cell[[marker_id]]$node_val, "numeric")) {
              max_abs_val <- max(abs(marker_limits[[marker_id]]))
              lims <- if (center_zero) c(-max_abs_val, max_abs_val) else c(0, max_abs_val)
              scale_color_gradientn(colours = colors, limits = lims)
            } else {
              scale_color_manual(values = colors)
            }
          }
        return(p)
      })

      if (flip) {
        p_markers <- wrap_plots(marker_plots, nrow = 1)
      } else {
        p_markers <- wrap_plots(marker_plots, ncol = 1)
      }

      return(p_markers)
    })

    # Wrap marker plots along rows or columns
    if (flip) {
      p <- wrap_plots(cell_plots, ncol = 1) & theme(legend.position = "top")
    } else {
      p <- wrap_plots(cell_plots, nrow = 1)
    }
    # Collect guides such that the color legend only appear on
    # the right side of the plot
    p <- p + plot_layout(guides = "collect")

    # row and col labels
    if (label_grid_axes && ((length(xyz_list_nested) > 1) || (length(xyz_list_nested[[1]]) > 1))) {
      # Create lists of row and column labels (grobs)
      # that can be added to the final patchwork
      marker_labels <- lapply(names(xyz_list_nested[[1]]), function(lbl) {
        wrap_elements(
          grid::textGrob(
            label = lbl
          )
        ) +
          theme(plot.background = element_rect(fill = "transparent", colour = NA))
      })
      cell_labels <- lapply(names(xyz_list_nested), function(lbl) {
        wrap_elements(
          grid::textGrob(
            label = lbl
          )
        ) +
          theme(plot.background = element_rect(fill = "transparent", colour = NA))
      })

      if (!flip) {
        col_labels <- cell_labels
        row_labels <- marker_labels
      } else {
        col_labels <- marker_labels
        row_labels <- cell_labels
      }

      col_width <- (1 - margin_widths[2]) / length(col_labels)
      row_width <- (1 - margin_widths[1]) / length(row_labels)

      # Create a patchwork with the column labels, leaving an empty
      # space for the upper left corner of the grid
      col_labels <- c(list(plot_spacer()), col_labels) %>%
        wrap_plots(
          nrow = 1,
          widths = c(margin_widths[2], rep(col_width, length(col_labels)))
        )

      row_labels <- wrap_plots(row_labels, ncol = 1)

      # Add row labels to the patchwork on the left side
      p <- wrap_plots(row_labels, p, widths = c(margin_widths[1], 1 - margin_widths[1]))
      # Add column labels to the patchwork on the top
      p <- wrap_plots(col_labels, p, ncol = 1, heights = c(margin_widths[2], 1 - margin_widths[2]))
    }
  }

  # Set background color to transparent so that the
  # background color for the frames can be controlled with
  # the png device
  p <- p & theme(plot.background = element_rect(fill = "transparent", color = NA))

  return(p)
}


#' Draws a scatter plot with simulated depth using base R
#'
#' @param xyz_list_nested A list of tibbles (\code{tbl_df}) with columns
#' 'x', 'y', 'z', and 'node_val'.
#' @param marker_limits A list with numeric vectors of length 2. The range
#' of the marker values.
#' @param colors A vector of valid colors.
#' @param center_zero A logical value indicating whether the plot should be
#' centered around zero.
#' @param xyz_limits A numeric vectors of length 2. The range of
#' the x, y, and z axes.
#' @param pt_opacity A numeric value between 0 and 1 indicating the opacity
#' of the points.
#' @param pt_size A numeric value indicating the maximum size of the points.
#' @param bg A valid color indicating the background color.
#' @param flip A logical value indicating whether the plot should be flipped.
#'
#' @noRd
#'
.scatter_3d_plot_base_R <- function(
  xyz_list_nested,
  marker_limits,
  colors,
  center_zero,
  xyz_limits,
  pt_opacity,
  pt_size,
  bg,
  flip = FALSE
) {
  # Save the current par settings
  # and reset when the function exits
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  if (flip) {
    ar <- c(length(xyz_list_nested[[1]]), length(xyz_list_nested))
    lyt <- matrix(1:(ar[1] * ar[2]), nrow = ar[2], byrow = FALSE)
  } else {
    ar <- c(length(xyz_list_nested), length(xyz_list_nested[[1]]))
    lyt <- matrix(1:(ar[1] * ar[2]), nrow = ar[1], byrow = TRUE)
  }

  # Create a grid for plotting
  par(
    mar = c(0, 0, 0, 0), # Remove margins
    pty = "s", # square aspect ratio <-> coord_fixed()
    bg = ifelse(bg == "transparent", NA, bg)
  )

  layout(mat = lyt)

  for (marker_id in names(xyz_list_nested[[1]])) {
    for (cell_id in names(xyz_list_nested)) {
      df <- xyz_list_nested[[cell_id]][[marker_id]]

      max_radius <- max(sqrt(rowSums((df %>% select(x, y, z))^2)))

      z_norm <- scales::rescale(df$z, from = c(-max_radius, max_radius), to = c(3, 1))
      apparent_sizes <- sqrt(1 / (z_norm^2))

      # Define node colors based on the type of node_val
      if (inherits(df$node_val, "numeric")) {
        if (center_zero) {
          max_abs_node_val <- max(marker_limits[[marker_id]])
          cols <- scales::col_numeric(domain = c(-max_abs_node_val, max_abs_node_val), palette = colors)(df$node_val)
        } else {
          cols <- scales::col_numeric(domain = marker_limits[[marker_id]], palette = colors)(df$node_val)
        }
      } else {
        if (inherits(df$node_val, "character")) {
          df$node_val <- as.factor(df$node_val)
        }
        cols <- scales::col_factor(domain = levels(df$node_val), palette = colors)(df$node_val)
      }

      # Add transparency to the colors
      cols <- scales::alpha(cols, alpha = pt_opacity)

      # Scale the apparent size of the points based on the y-axis values
      # plot() uses the cex parameter to scale the size of the points
      # which is based on the radius of the circle.
      # First, we need to normalize the y-axis values to a range of 3 to 1,
      # meaning that the point closest to the camera is 1 distance away and the
      # point furthest away is 3 distances away.
      y_norm <- scales::rescale(df$y, from = c(-max_radius, max_radius), to = c(3, 1))
      # Then we use the inverse square law to adjust the area of the points
      # and take the square root to get the radii. By multiplying the radii
      # by pt_size, we can adjust the size of the points while retaining the
      # perceived size differences.
      apparent_sizes <- sqrt(1 / (y_norm^2)) * pt_size

      # TODO: Add column and row labels
      plot(
        df$x,
        df$z,
        cex = apparent_sizes,
        pch = 16,
        col = cols,
        xlab = "",
        ylab = "",
        axes = FALSE,
        xlim = xyz_limits,
        ylim = xyz_limits,
      )
    }
  }
}

#' A helper function to create a rotation matrix
#' for a given angle in degrees. The rotation is
#' around the z-axis.
#'
#' @param angle_degrees A numeric value indicating the angle of rotation
#' in degrees.
#'
#' @returns A 3x3 matrix representing the rotation matrix.
#'
#' @noRd
#'
.create_rotation_matrix <- function(angle_degrees) {
  angle_radians <- angle_degrees * pi / 180
  matrix(
    c(
      cos(angle_radians), -sin(angle_radians), 0,
      sin(angle_radians), cos(angle_radians), 0,
      0, 0, 1
    ),
    nrow = 3, byrow = TRUE
  )
}


#' A helper function to rotate a set of 3D coordinates
#' around the z-axis.
#'
#' @param coords A tibble with the x, y, and z coordinates
#' @param angle_degrees A numeric value indicating the angle of rotation
#' in degrees.
#'
#' @returns A tibble with the rotated coordinates.
#'
#' @noRd
#'
.rotate_around_z <- function(coords, angle_degrees) {
  rotation_matrix <- .create_rotation_matrix(angle_degrees)
  as_tibble(t(rotation_matrix %*% t(coords %>% as.matrix())), .name_repair = ~ c("x", "y", "z"))
}


#' Scale layout coordinates
#'
#' This function scales the layout coordinates to a given percentile
#' of the maximum radius.
#'
#' @param layout A tibble with x, y and optionally z coordinates.
#' @param percentile A numeric value between 0 and 1 indicating the
#' percentile of the maximum radius to scale the layout to.
#'
#' @return A tibble with the scaled coordinates.
#'
#' @export
#'
scale_layout <- function(
  layout,
  percentile = 0.95
) {
  assert_class(layout, classes = c("data.frame", "tbl_df"))
  assert_within_limits(percentile, c(0, 1))

  layout <- layout %>%
    select(any_of(c("x", "y", "z")))

  if (!ncol(layout) %in% c(2, 3)) {
    abort("The layout must have 2 or 3 columns.")
  }

  radii <- sqrt(rowSums(layout^2))
  perc <- quantile(radii, percentile)
  layout_scaled <- layout %>%
    mutate(across(everything(), ~ .x / perc))

  return(layout_scaled)
}

#' A helper function to render the first frame of the video
#'
#' In interactive sessions, the first frame of the video is rendered
#' in the viewer. The user can then decide whether to continue
#' rendering the rest of the frames.
#'
#' In non-interactive sessions, the path to the first frame is
#' returned. The user can then decide whether to continue
#' rendering the rest of the frames.
#'
#' @param xyz_list_nested A list of tibbles with the data to plot
#' @param marker_limits A numeric vector of length 2 with the limits
#' @param pt_opacity A numeric value with the opacity of the points
#' @param pt_size A numeric value with the size of the points
#' @param marker_col,cell_col A character string with the column to facet by
#' @param xyz_limits A numeric vector of length 2 with the range
#' of the x, y, and z axes
#' @param colors A character vector with the colors to use
#' @param center_zero A logical indicating whether to center the
#' color scale at zero
#' @param pad A numeric value with the amount of padding to add
#' @param tmp_dir A character string with the path to the temporary
#' directory
#' @param width,height A numeric value specifying the width and
#' height of the plot
#' @param res A numeric value specifying the resolution of the plot
#' @param ggplot_theme A ggplot2 theme object
#' @param title A character string with the title of the plot
#' @param bg A character string with the background color of the plot
#' @param use_facet_grid A logical indicating whether to use facet_grid
#' @param flip A logical indicating whether to flip the plot
#' @param label_grid_axes A logical indicating whether to label the axes
#' @param margin_widths A numeric vector of length 2 with the margin widths
#' @param graphics_use A character string indicating the graphics library
#' to use
#' @param dev_png A valid png device, e.g. \code{ragg::agg_png()}
#'
#' @returns Nothing
#'
#' @noRd
#'
.render_first_frame <- function(
  xyz_list_nested,
  marker_limits,
  pt_opacity,
  pt_size,
  marker_col,
  cell_col,
  colors,
  center_zero,
  xyz_limits,
  tmp_dir,
  width,
  height,
  res,
  ggplot_theme,
  title,
  bg,
  use_facet_grid,
  flip,
  label_grid_axes,
  margin_widths,
  graphics_use,
  dev_png
) {
  rlang::check_installed("png")
  file <- file.path(tmp_dir, "plot_0000.png")

  # Use the more flexible but slower ggplot2 graphics
  if (graphics_use == "ggplot2") {
    p <- .scatter_3d_plot(
      xyz_list_nested,
      cell_col,
      marker_col,
      marker_limits,
      colors,
      center_zero,
      xyz_limits,
      pt_opacity,
      pt_size,
      use_facet_grid,
      flip,
      label_grid_axes,
      margin_widths
    )

    # Add theme
    if (!is.null(ggplot_theme)) {
      p <- p & ggplot_theme
    }
    if (title != "") {
      p <- p + plot_annotation(title = title)
    }

    # Export first frame
    ggsave(file, p, width = width, height = height, dpi = res, device = dev_png, units = "px", bg = bg)
  }

  # Use the faster but less flexible base graphics
  if (graphics_use == "base") {
    dev_png(file, width = width, height = height, res = res, units = "px", bg = bg)
    .scatter_3d_plot_base_R(
      xyz_list_nested,
      marker_limits,
      colors,
      center_zero,
      xyz_limits,
      pt_opacity,
      pt_size,
      bg,
      flip
    )
    dev.off()
  }

  # Display the first frame in the viewer
  p <- png::readPNG(file) %>% as.raster()
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mar = c(0, 0, 0, 0))
  plot(p)
  cli_alert(glue(
    "The first frame has been rendered in ",
    "the viewer. Do you want to continue?"
  ))

  # Ask user if they want to continue
  if (utils::menu(c("Yes", "No")) != 1) {
    abort("Aborted by user", call = caller_env())
  }

  dev.off()
  fs::file_delete(file)
}


#' Utility function to validate colors
#'
#' @param colors A vector of color named or hex values
#'
#' @returns A logical vector
#'
#' @examples
#' # Valid colors
#' .areColors(c("red", "steelblue", "#000", "#808080", "#80808080"))
#'
#' # Invalid colors
#' .areColors(c("RED"))
#'
#' @noRd
#'
.areColors <- function(colors) {
  sapply(colors, function(X) {
    tryCatch(is.matrix(col2rgb(X)),
             error = function(e) FALSE
    )
  })
}


#' Utility function to validate input parameters for render_rotating_layout
#'
#' @param data A tibble object
#' @param marker_col,cell_col A character string
#' @param pt_opacity A numeric value between 0 and 1
#' @param pt_size A numeric value
#' @param colors A vector of color names
#' @param max_degree A numeric value between 90 and 360
#' @param center_zero A logical value
#' @param scale_layout A logical value
#' @param frames An integer value
#' @param pad A numeric value
#' @param show_first_frame A logical value
#' @param width,height An integer value
#' @param res A numeric value
#' @param delay A numeric value
#' @param ggplot_theme A ggplot2 theme object
#' @param title A character string
#' @param bg A character string
#' @param label_grid_axes A logical value
#' @param margin_widths A numeric vector
#' @param use_facet_grid A logical value
#' @param flip A logical value
#' @param boomerang A logical value
#'
#' @returns Nothing
#'
#' @noRd
#'
.validate_render_rotating_layout_input_params <- function(
  data,
  marker_col,
  cell_col,
  pt_opacity,
  pt_size,
  colors,
  max_degree,
  center_zero,
  scale_layout,
  frames,
  pad,
  show_first_frame,
  width,
  height,
  res,
  delay,
  ggplot_theme,
  title,
  bg,
  label_grid_axes,
  margin_widths,
  use_facet_grid,
  flip,
  boomerang,
  call = caller_env()
) {
  assert_class(data, "tbl_df", call = call)
  assert_within_limits(pt_opacity, c(0, 1), call = call)
  assert_within_limits(pt_size, c(0, 5), call = call)
  assert_vector(colors, "character", call = call)
  assert_within_limits(max_degree, c(90, 360), call = call)
  assert_single_value(center_zero, "bool", call = call)
  assert_single_value(scale_layout, "bool", call = call)
  assert_single_value(frames, "integer", call = call)
  assert_within_limits(pad, c(0, 1), call = call)
  assert_single_value(show_first_frame, "bool", call = call)
  assert_single_value(width, "integer", call = call)
  assert_single_value(height, "integer", call = call)
  assert_single_value(res, "integer", call = call)
  assert_within_limits(delay, c(0.01, 1), call = call)
  assert_single_value(title, "string", call = call)
  assert_single_value(label_grid_axes, "bool", call = call)
  assert_length(margin_widths, n = 2, call = call)
  assert_within_limits(margin_widths, c(0, 1), call = call)
  assert_single_value(use_facet_grid, "bool", call = call)
  assert_single_value(flip, "bool", call = call)
  assert_single_value(boomerang, "bool", call = call)

  if (!all(.areColors(colors))) {
    cli::cli_abort(
      c(
        "x" = "Invalid color values in {.var colors}.",
        "i" = "All colors must be valid color strings."
      ),
      call = call
    )
  }
  if (!(length(bg) == 1 && .areColors(bg))) {
    cli::cli_abort(
      c(
        "x" = "{.var bg} must be a valid color."
      ),
      call = call
    )
  }

  missing_cols <- setdiff(c("x", "y", "z", "node_val"), names(data))
  if (length(missing_cols) > 0) {
    cli::cli_abort(
      c(
        "x" = "Missing column{?s}: {.str {missing_cols}}"
      ),
      call = call
    )
  }

  # Check marker_col and cell_col
  if (!is.null(marker_col)) {
    assert_single_value(marker_col, "string", call = call)
    assert_col_in_data(marker_col, data, call = call)
  }
  if (!is.null(cell_col)) {
    assert_single_value(cell_col, "string", call = call)
    assert_col_in_data(cell_col, data, call = call)
  }
  if (!is.null(marker_col) && !is.null(cell_col)) {
    assert_single_values_are_different(marker_col, cell_col, call = call)
  }

  # Check column types
  assert_col_class("x", data, "numeric", call = call)
  assert_col_class("y", data, "numeric", call = call)
  assert_col_class("z", data, "numeric", call = call)

  if (inherits(data$node_val, c("character", "factor"))) {
    if (length(colors) < length(unique(data$node_val))) {
      cli::cli_warn("Not enough colors for node_val groups. Expanding colors vector.")
      colors <- colorRampPalette(colors)(length(unique(data$node_val)))
    }
    if (inherits(data$node_val, "character")) {
      data <- data %>%
        mutate(node_val = factor(node_val, levels = unique(node_val)))
    }
  }
  if (inherits(data$node_val, "numeric")) {
    if (any(is.na(data$node_val))) {
      cli::cli_abort(
        c(
          "x" = "Column {.val node_val} cannot contain missing values"
        ),
        call = call
      )
    }
  }
  if (!is.null(ggplot_theme)) {
    assert_class(ggplot_theme, "theme", call = call)
  }
}
