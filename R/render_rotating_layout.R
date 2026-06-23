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
#' @section illumination:
#' When \code{use_illumination = TRUE}, node colors are derived from \code{node_val}
#' and then modulated by a geometry-based illumination mask computed with
#' \code{\link{heuristic_illumination}}. This adds depth cues on top of the marker
#' color encoding. By default, shadows are blended in HSV space using
#' \code{illumination_ambient} as the minimum brightness floor (not the
#' ambient-occlusion weight in \code{heuristic_illumination}) and
#' \code{illumination_sat_boost} for saturation compensation in shadowed regions.
#' When \code{illumination_shadow_colors} is set, shadows are instead blended
#' toward a second color gradient mapped from the illumination mask; in that mode
#' \code{illumination_sat_boost} is ignored. \code{PixelgenGradient(n, "NaturalBlue")}
#' works well as a shadow palette while \code{colors} carries the marker signal.
#' Set \code{normalize_illumination = FALSE} to use raw output from
#' \code{heuristic_illumination} instead of rescaling the mask to \code{[0, 1]}.
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
#' @param max_degree A numeric value between 0 and 360. The maximum angle
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
#' @param use_illumination A logical value indicating whether to modulate
#' \code{node_val} colors with a heuristic illumination mask. Default is \code{FALSE}.
#' @param illumination_ambient A numeric value between 0 and 1 specifying the
#' minimum brightness floor when blending illumination. Default is \code{0.3}.
#' @param illumination_sat_boost A non-negative numeric value controlling
#' saturation compensation in shadowed regions when using HSV blending.
#' Ignored when \code{illumination_shadow_colors} is set. Increase this to increase
#' saturation in darker points. Default is \code{0.6}.
#' @param illumination_shadow_colors \code{NULL} or a vector of valid colors.
#' When \code{NULL} (default), illumination is blended in HSV space. When set,
#' each point's illumination is mapped to this palette instead and linearly interpolated
#' with the base \code{node_val} color.
#' @param normalize_illumination A logical value indicating whether the
#' illumination mask should be rescaled to \code{[0, 1]} per cell. Default is
#' \code{TRUE}.
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
  keep_frames = FALSE,
  use_illumination = FALSE,
  illumination_ambient = 0.3,
  illumination_sat_boost = 0.6,
  illumination_shadow_colors = NULL,
  normalize_illumination = TRUE
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
    flip, boomerang, use_illumination, illumination_ambient,
    illumination_sat_boost, illumination_shadow_colors, normalize_illumination
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

  if (use_illumination) {
    # Compute shading once per cell.
    xyz_list <- lapply(xyz_list, function(xyz) {
      ill <- heuristic_illumination(xyz)
      if (normalize_illumination) {
        ill <- scales::rescale(ill, to = c(0, 1))
      }
      xyz$illumination <- ill
      xyz
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

  if (use_illumination) {
    # Precompute colors before the frame loop; only coordinates rotate per frame.
    xyz_list_nested <- lapply(xyz_list_nested, function(xyz_list_cell) {
      lapply(names(xyz_list_cell), function(marker_id) {
        df <- xyz_list_cell[[marker_id]]
        df$point_color <- .compute_illuminated_point_colors(
          df,
          colors,
          marker_limits,
          marker_id,
          center_zero,
          illumination_ambient,
          illumination_sat_boost,
          illumination_shadow_colors,
          normalize_illumination
        )
        df
      }) %>%
        set_names(names(xyz_list_cell))
    }) %>%
      set_names(names(xyz_list_nested))
  }

  # Show the first frame if requested
  if (show_first_frame && interactive()) {
    .render_first_frame(
      xyz_list_nested, marker_limits, pt_opacity, pt_size,
      marker_col, cell_col, colors, center_zero, xyz_limits,
      tmp_dir, width, height, res, ggplot_theme, title, bg,
      use_facet_grid, flip, label_grid_axes, margin_widths,
      graphics_use, dev_png, use_illumination
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
          margin_widths,
          use_illumination
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
          flip,
          use_illumination
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


#' Blend base hex colors with an illumination mask in HSV space
#'
#' @param hex_colors A character vector of hex colors.
#' @param illum_mask A numeric vector of illumination values in `[0, 1]`.
#' @param ambient_intensity Minimum brightness floor in `[0, 1]`.
#' @param sat_boost Non-negative saturation boost in shadowed regions.
#'
#' @returns A character vector of hex colors.
#'
#' @noRd
#'
.apply_hsv_illumination <- function(
  hex_colors,
  illum_mask,
  ambient_intensity = 0.1,
  sat_boost = 0.7,
  call = caller_env()
) {
  assert_valid_color(hex_colors, n = 1, call = call)
  assert_vector(illum_mask, type = "numeric", n = 1, call = call)
  assert_vectors_x_y_length_equal(hex_colors, illum_mask, call = call)
  assert_within_limits(illum_mask, c(0, 1), call = call)
  assert_single_value(ambient_intensity, type = "numeric", call = call)
  assert_within_limits(ambient_intensity, c(0, 1), call = call)
  assert_single_value(sat_boost, type = "numeric", call = call)
  assert_within_limits(sat_boost, c(0, Inf), call = call)

  rgb_mat <- grDevices::col2rgb(hex_colors)
  hsv_mat <- grDevices::rgb2hsv(rgb_mat)

  h <- hsv_mat[1, ]
  s <- hsv_mat[2, ]
  v <- hsv_mat[3, ]

  light_factor <- ambient_intensity + (illum_mask * (1 - ambient_intensity))
  v_new <- v * light_factor
  s_boosted <- s * (1 + (1 - light_factor) * sat_boost)
  s_new <- pmin(1, s_boosted)

  grDevices::hsv(h = h, s = s_new, v = v_new)
}


#' Blend base hex colors toward shadow colors using RGB interpolation
#'
#' @param base_hex A character vector of base hex colors.
#' @param shadow_hex A character vector of shadow hex colors (length 1 or
#'   equal to \code{base_hex}).
#' @param illum_mask A numeric vector of illumination values in `[0, 1]`.
#' @param ambient_intensity Minimum light factor in `[0, 1]`.
#'
#' @returns A character vector of hex colors.
#'
#' @noRd
#'
.apply_palette_illumination <- function(
  base_hex,
  shadow_hex,
  illum_mask,
  ambient_intensity = 0.1,
  call = caller_env()
) {
  assert_valid_color(base_hex, n = 1, call = call)
  assert_valid_color(shadow_hex, n = 1, call = call)
  assert_vector(illum_mask, type = "numeric", n = 1, call = call)
  assert_vectors_x_y_length_equal(base_hex, illum_mask, call = call)
  assert_within_limits(illum_mask, c(0, 1), call = call)
  if (length(shadow_hex) != 1L) {
    assert_vectors_x_y_length_equal(base_hex, shadow_hex, call = call)
  }
  assert_single_value(ambient_intensity, type = "numeric", call = call)
  assert_within_limits(ambient_intensity, c(0, 1), call = call)

  rgb_base <- grDevices::col2rgb(base_hex)
  if (length(shadow_hex) == 1L) {
    rgb_shadow <- matrix(grDevices::col2rgb(shadow_hex), nrow = 3, ncol = length(base_hex))
  } else {
    rgb_shadow <- grDevices::col2rgb(shadow_hex)
  }

  light_factor <- ambient_intensity + (illum_mask * (1 - ambient_intensity))
  mask_mat <- matrix(light_factor, nrow = 3, ncol = length(light_factor), byrow = TRUE)
  inv_mask_mat <- matrix(1 - light_factor, nrow = 3, ncol = length(light_factor), byrow = TRUE)

  rgb_blended <- (rgb_base * mask_mat) + (rgb_shadow * inv_mask_mat)

  grDevices::rgb(
    red = rgb_blended[1, ],
    green = rgb_blended[2, ],
    blue = rgb_blended[3, ],
    maxColorValue = 255
  )
}


#' Compute illuminated point colors from node values and illumination
#'
#' @param df A tibble with \code{node_val} and \code{illumination} columns.
#' @param colors A vector of valid colors.
#' @param marker_limits A list with numeric vectors of length 2, or \code{NULL}.
#' @param marker_id A character string with the marker identifier.
#' @param center_zero A logical value indicating whether to center the scale at zero.
#' @param illumination_ambient Minimum brightness floor in `[0, 1]`.
#' @param illumination_sat_boost Non-negative saturation boost in shadowed regions.
#' @param illumination_shadow_colors \code{NULL} or a vector of shadow palette colors.
#' @param normalize_illumination A logical value indicating whether to rescale
#' illumination to \code{[0, 1]}.
#'
#' @returns A character vector of hex colors.
#'
#' @noRd
#'
.compute_illuminated_point_colors <- function(
  df,
  colors,
  marker_limits,
  marker_id,
  center_zero,
  illumination_ambient,
  illumination_sat_boost,
  illumination_shadow_colors = NULL,
  normalize_illumination = TRUE,
  call = caller_env()
) {
  assert_class(df, c("data.frame", "tbl_df"), call = call)
  assert_col_in_data("node_val", df, call = call)
  assert_col_in_data("illumination", df, call = call)
  assert_col_class("illumination", df, "numeric", call = call)
  assert_single_value(marker_id, type = "string", call = call)
  assert_single_value(center_zero, type = "bool", call = call)
  assert_single_value(illumination_ambient, "numeric", call = call)
  assert_within_limits(illumination_ambient, c(0, 1), call = call)
  assert_single_value(illumination_sat_boost, "numeric", call = call)
  assert_within_limits(illumination_sat_boost, c(0, Inf), call = call)
  assert_vector(colors, type = "character", n = 1, call = call)
  assert_valid_color(colors, n = 1, call = call)
  assert_valid_color(illumination_shadow_colors, n = 1, allow_null = TRUE, call = call)
  assert_single_value(normalize_illumination, "bool", call = call)

  if (inherits(df$node_val, "numeric")) {
    assert_class(marker_limits, "list", call = call)
    assert_x_in_y(marker_id, names(marker_limits), call = call)
    lims <- .node_val_lims(df$node_val, marker_limits, marker_id, center_zero)
    base_color <- scales::col_numeric(
      palette = colors,
      domain = lims,
      na.color = "transparent"
    )(df$node_val)
  } else {
    node_val <- df$node_val
    if (inherits(node_val, "character")) {
      node_val <- factor(node_val, levels = unique(node_val))
    }

    if (!is.null(names(colors))) {
      assert_x_in_y(levels(node_val), names(colors), call = call)

      colors <-
        colors[match(levels(node_val), names(colors))]
    }

    base_color <- scales::col_factor(
      domain = levels(node_val),
      palette = colors
    )(node_val)
  }

  illum_mask <- df$illumination

  if (is.null(illumination_shadow_colors)) {
    return(.apply_hsv_illumination(
      base_color,
      illum_mask,
      ambient_intensity = illumination_ambient,
      sat_boost = illumination_sat_boost,
      call = call
    ))
  }

  shadow_color <- scales::col_numeric(
    palette = illumination_shadow_colors,
    domain = if (normalize_illumination) {
      c(0, 1)
    } else {
      range(illum_mask, na.rm = TRUE)
    },
    na.color = "transparent"
  )(illum_mask)

  .apply_palette_illumination(
    base_color,
    shadow_color,
    illum_mask,
    ambient_intensity = illumination_ambient,
    call = call
  )
}


#' Numeric limits for node_val color scales
#'
#' @param node_val A numeric vector of node values.
#' @param marker_limits A list with numeric min/max vectors, or \code{NULL}.
#' @param marker_id A character string with the marker identifier, or \code{NULL}.
#' @param center_zero A logical value indicating whether to center the scale at zero.
#'
#' @returns A numeric vector of length 2.
#'
#' @noRd
#'
.node_val_lims <- function(
  node_val,
  marker_limits = NULL,
  marker_id = NULL,
  center_zero = FALSE
) {
  if (!is.null(marker_id) && !is.null(marker_limits)) {
    lims <- marker_limits[[marker_id]]
  } else {
    lims <- range(node_val, na.rm = TRUE)
  }
  if (center_zero) {
    max_abs_val <- max(abs(lims))
    lims <- c(-max_abs_val, max_abs_val)
  }
  lims
}


#' Minimal tibble to drive the node_val color legend when colors are precomputed
#'
#' @param node_val A numeric or factor vector of node values.
#' @param marker_limits A list with numeric vectors of length 2, or \code{NULL}.
#' @param marker_id A character string with the marker identifier, or \code{NULL}.
#' @param center_zero A logical value indicating whether to center the scale at zero.
#'
#' @returns A tibble with columns \code{node_val}, \code{x}, \code{y}, and \code{z}.
#'
#' @noRd
#'
.legend_carrier_data <- function(
  node_val,
  marker_limits = NULL,
  marker_id = NULL,
  center_zero = FALSE
) {
  if (inherits(node_val, "numeric")) {
    lims <- .node_val_lims(node_val, marker_limits, marker_id, center_zero)
    tibble::tibble(
      node_val = lims,
      x = 0,
      y = 0,
      z = 0
    )
  } else {
    if (inherits(node_val, "character")) {
      node_val <- factor(node_val, levels = unique(node_val))
    }
    tibble::tibble(
      node_val = factor(levels(node_val), levels = levels(node_val)),
      x = 0,
      y = 0,
      z = 0
    )
  }
}


#' Create a ggplot2 color scale for node values
#'
#' @param node_val A numeric or factor vector of node values.
#' @param colors A vector of valid colors.
#' @param center_zero A logical value indicating whether to center the scale at zero.
#' @param marker_limits A list with numeric vectors of length 2, or \code{NULL}.
#' @param marker_id A character string with the marker identifier, or \code{NULL}.
#' @param use_illumination A logical value indicating whether point colors are
#' precomputed outside the color aesthetic.
#' @param pt_size Maximum point size used for legend key sizing.
#'
#' @returns A ggplot2 color scale.
#'
#' @noRd
#'
.node_val_color_scale <- function(
  node_val,
  colors,
  center_zero,
  marker_limits = NULL,
  marker_id = NULL,
  use_illumination = FALSE,
  pt_size = 1
) {
  if (inherits(node_val, "numeric")) {
    lims <- .node_val_lims(node_val, marker_limits, marker_id, center_zero)
    scale_color_gradientn(colours = colors, limits = lims)
  } else {
    if (inherits(node_val, "character")) {
      node_val <- factor(node_val, levels = unique(node_val))
    }
    # Colors are precomputed outside aes; override.aes keeps legend keys visible.
    # "legend" is ggplot2's default guide shorthand (same as omitting guide=).
    color_guide <- if (use_illumination) {
      ggplot2::guide_legend(override.aes = list(alpha = 1, size = pt_size))
    } else {
      "legend"
    }
    scale_color_manual(
      values = colors,
      limits = levels(node_val),
      guide = color_guide
    )
  }
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
#' @param use_illumination A logical value indicating whether precomputed
#' \code{point_color} values should be used.
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
  margin_widths = c(0.1, 0.1),
  use_illumination = FALSE
) {
  if (use_facet_grid) {
    xyz_collapsed <- lapply(xyz_list_nested, function(xyz_cell) {
      xyz_cell %>%
        unname() %>%
        bind_rows()
    }) %>%
      unname() %>%
      bind_rows()
    if (use_illumination) {
      legend_df <- .legend_carrier_data(
        xyz_collapsed$node_val,
        marker_limits,
        center_zero = center_zero
      )
      p <- ggplot(xyz_collapsed, aes(x, z, size = y)) +
        geom_point(color = xyz_collapsed$point_color, alpha = pt_opacity) +
        geom_point(
          data = legend_df,
          aes(x, z, color = node_val),
          inherit.aes = FALSE,
          alpha = 0,
          size = 0
        ) +
        scale_size(range = c(0, pt_size), limits = xyz_limits) +
        coord_fixed() +
        theme_void() +
        guides(size = "none") +
        scale_x_continuous(limits = xyz_limits, expand = expansion()) +
        scale_y_continuous(limits = xyz_limits, expand = expansion()) +
        .node_val_color_scale(
          xyz_collapsed$node_val,
          colors,
          center_zero,
          marker_limits = marker_limits,
          use_illumination = use_illumination,
          pt_size = pt_size
        ) +
        {
          if (flip) {
            facet_grid(reformulate(marker_col, cell_col), switch = "y")
          } else {
            facet_grid(reformulate(cell_col, marker_col), switch = "y")
          }
        }
    } else {
      p <- ggplot(xyz_collapsed, aes(x, z, size = y, color = node_val)) +
        geom_point(alpha = pt_opacity) +
        scale_size(range = c(0, pt_size), limits = xyz_limits) +
        coord_fixed() +
        theme_void() +
        guides(size = "none") +
        scale_x_continuous(limits = xyz_limits, expand = expansion()) +
        scale_y_continuous(limits = xyz_limits, expand = expansion()) +
        .node_val_color_scale(
          xyz_collapsed$node_val,
          colors,
          center_zero,
          marker_limits = marker_limits,
          use_illumination = use_illumination,
          pt_size = pt_size
        ) +
        {
          if (flip) {
            facet_grid(reformulate(marker_col, cell_col), switch = "y")
          } else {
            facet_grid(reformulate(cell_col, marker_col), switch = "y")
          }
        }
    }
  }

  if (!use_facet_grid) {
    # Create a list of ggplot objects per cell
    cell_plots <- lapply(names(xyz_list_nested), function(cell_id) {
      xyz_list_cell <- xyz_list_nested[[cell_id]]
      # Create a list of ggplot objects per marker
      marker_plots <- lapply(names(xyz_list_cell), function(marker_id) {
        df <- xyz_list_cell[[marker_id]] %>%
          mutate(row_text = marker_id)
        if (use_illumination) {
          legend_df <- .legend_carrier_data(
            df$node_val,
            marker_limits,
            marker_id,
            center_zero
          )
          p <- ggplot(df, aes(x, z, size = y)) +
            geom_point(color = df$point_color, alpha = pt_opacity) +
            geom_point(
              data = legend_df,
              aes(x, z, color = node_val),
              inherit.aes = FALSE,
              alpha = 0,
              size = 0
            ) +
            scale_size(range = c(0, pt_size), limits = xyz_limits) +
            coord_fixed() +
            theme_void() +
            guides(size = "none") +
            scale_x_continuous(limits = xyz_limits, expand = expansion()) +
            scale_y_continuous(limits = xyz_limits, expand = expansion()) +
            .node_val_color_scale(
              df$node_val,
              colors,
              center_zero,
              marker_limits = marker_limits,
              marker_id = marker_id,
              use_illumination = use_illumination,
              pt_size = pt_size
            )
        } else {
          p <- ggplot(df, aes(x, z, size = y, color = node_val)) +
            geom_point(alpha = pt_opacity) +
            scale_size(range = c(0, pt_size), limits = xyz_limits) +
            coord_fixed() +
            theme_void() +
            guides(size = "none") +
            scale_x_continuous(limits = xyz_limits, expand = expansion()) +
            scale_y_continuous(limits = xyz_limits, expand = expansion()) +
            .node_val_color_scale(
              df$node_val,
              colors,
              center_zero,
              marker_limits = marker_limits,
              marker_id = marker_id,
              use_illumination = use_illumination,
              pt_size = pt_size
            )
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
#' @param use_illumination A logical value indicating whether precomputed
#' \code{point_color} values should be used.
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
  flip = FALSE,
  use_illumination = FALSE
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

  graphics::layout(mat = lyt)

  for (marker_id in names(xyz_list_nested[[1]])) {
    for (cell_id in names(xyz_list_nested)) {
      df <- xyz_list_nested[[cell_id]][[marker_id]]

      max_radius <- max(sqrt(rowSums((df %>% select(x, y, z))^2)))

      z_norm <- scales::rescale(df$z, from = c(-max_radius, max_radius), to = c(3, 1))
      apparent_sizes <- sqrt(1 / (z_norm^2))

      # Define node colors based on the type of node_val
      if (use_illumination) {
        cols <- df$point_color
      } else if (inherits(df$node_val, "numeric")) {
        if (center_zero) {
          max_abs_node_val <- max(abs(marker_limits[[marker_id]]))
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
#' @param use_illumination A logical value indicating whether precomputed
#' \code{point_color} values should be used.
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
  dev_png,
  use_illumination = FALSE
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
      margin_widths,
      use_illumination
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
      flip,
      use_illumination
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
#' @param max_degree A numeric value between 0 and 360
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
#' @param use_illumination A logical value
#' @param illumination_ambient A numeric value
#' @param illumination_sat_boost A numeric value
#' @param illumination_shadow_colors A character vector or \code{NULL}
#' @param normalize_illumination A logical value
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
  use_illumination,
  illumination_ambient,
  illumination_sat_boost,
  illumination_shadow_colors,
  normalize_illumination,
  call = caller_env()
) {
  assert_class(data, "tbl_df", call = call)
  assert_within_limits(pt_opacity, c(0, 1), call = call)
  assert_within_limits(pt_size, c(0, 5), call = call)
  assert_vector(colors, "character", call = call)
  assert_within_limits(max_degree, c(0, 360), call = call)
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
  assert_single_value(use_illumination, "bool", call = call)
  assert_single_value(normalize_illumination, "bool", call = call)
  assert_single_value(illumination_ambient, "numeric", call = call)
  assert_within_limits(illumination_ambient, c(0, 1), call = call)
  assert_single_value(illumination_sat_boost, "numeric", call = call)
  assert_within_limits(illumination_sat_boost, c(0, Inf), call = call)

  if (isTRUE(use_illumination) && !is.null(illumination_shadow_colors)) {
    assert_valid_color(illumination_shadow_colors, n = 1, call = call)
  }

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

#' Compute heuristic illumination for 3D layouts
#'
#' Combines three simple lighting heuristics for 3D coordinates:
#' (1) directional light from the positive z-axis,
#' (2) radial volume shading from the origin,
#' and (3) ambient occlusion approximated from mean distance to nearest neighbors.
#'
#' @param layout A data frame or tibble with numeric columns `x`, `y`, and `z`.
#' @param clamp_quantiles Numeric vector of length 2 in `[0, 1]`. Illumination is
#'   clamped to these quantiles to reduce outlier influence. Default: `c(0.01, 0.95)`.
#' @param directional_light_weight Non-negative numeric scalar. Weight for directional
#'   light component. Default: `0.7`.
#' @param volume_shading_weight Non-negative numeric scalar. Weight for radial volume
#'   shading component. Default: `0.5`.
#' @param ambient_occlusion_weight Non-negative numeric scalar. Weight for ambient
#'   occlusion component. Default: `1`.
#' @param ambient_occlusion_k Positive integer. Number of nearest neighbors used for
#'   ambient occlusion approximation. Default: `20`.
#' @param normalize_weights Logical; if `TRUE`, weights are normalized to sum to 1.
#'   Default: `TRUE`.
#'
#' @returns A numeric vector of illumination values (length `nrow(layout)`). Higher values indicate stronger
#'   illumination.
#'
#' @examples
#'
#' library(dplyr)
#' set.seed(1)
#'
#' # Here we simulate some 3D coordinates with a roughly spherical distribution
#' n <- 20000
#' n_surface <- 19000
#' n_interior <- 1000
#'
#' # Surface points: normalize to unit sphere, add small Gaussian noise
#' xyz_surface <- matrix(rnorm(n_surface * 3), ncol = 3)
#' xyz_surface <- xyz_surface / sqrt(rowSums(xyz_surface^2)) # project to unit sphere
#' xyz_surface <- xyz_surface + matrix(rnorm(n_surface * 3, sd = 0.05), ncol = 3)
#'
#' # Interior points: uniform in ball via rejection sampling
#' xyz_interior <- matrix(rnorm(n_interior * 3), ncol = 3)
#' radii <- runif(n_interior)^(1 / 3) # cube root for uniform volume distribution
#' xyz_interior <- xyz_interior / sqrt(rowSums(xyz_interior^2)) * radii * 0.8
#'
#' layout <- tibble::tibble(
#'   x = c(xyz_surface[, 1], xyz_interior[, 1]),
#'   y = c(xyz_surface[, 2], xyz_interior[, 2]),
#'   z = c(xyz_surface[, 3], xyz_interior[, 3])
#' )
#' illum <- heuristic_illumination(layout)
#'
#' # Create a temporary GIF file and render a rotating layout
#' # using the computed illumination as node values
#' temp_gif <- fs::file_temp(ext = ".gif")
#' render_rotating_layout(
#'   data = layout %>%
#'     mutate(node_val = illum),
#'   pt_size = 0.8,
#'   width = 740,
#'   height = 650,
#'   colors = PixelgenGradient(100, "NaturalBlue"),
#'   file = temp_gif,
#'   max_degree = 30,
#'   frames = 20,
#'   delay = 1 / 20,
#'   res = 100,
#'   boomerang = TRUE,
#'   show_first_frame = FALSE
#' )
#'
#' @export
heuristic_illumination <- function(
  layout,
  clamp_quantiles = c(0.01, 0.95),
  directional_light_weight = 0.7,
  volume_shading_weight = 0.5,
  ambient_occlusion_weight = 1,
  ambient_occlusion_k = 20,
  normalize_weights = TRUE
) {
  expect_FNN()

  assert_class(layout, c("data.frame", "tbl_df"))
  assert_x_in_y(x = c("x", "y", "z"), y = names(layout))
  assert_vector(layout$x, "numeric")
  assert_vector(layout$y, "numeric")
  assert_vector(layout$z, "numeric")

  coords <- as.matrix(layout[, c("x", "y", "z")])

  if (any(!is.finite(coords))) {
    cli::cli_abort("Columns `x`, `y`, and `z` must contain only finite values.")
  }

  assert_vector(clamp_quantiles, "numeric", n = 2)
  assert_within_limits(clamp_quantiles, c(0, 1))

  if (clamp_quantiles[1] >= clamp_quantiles[2]) {
    cli::cli_abort("`clamp_quantiles[1]` must be less than `clamp_quantiles[2]`.")
  }

  assert_single_value(directional_light_weight, "numeric")
  assert_within_limits(directional_light_weight, c(0, Inf))
  assert_single_value(volume_shading_weight, "numeric")
  assert_within_limits(volume_shading_weight, c(0, Inf))
  assert_single_value(ambient_occlusion_weight, "numeric")
  assert_within_limits(ambient_occlusion_weight, c(0, Inf))
  assert_single_value(ambient_occlusion_k, "integer")
  assert_within_limits(ambient_occlusion_k, c(1, nrow(layout) - 1))
  assert_single_value(normalize_weights, "bool")

  # Rescale function
  safe_rescale <- function(x, to = c(0, 1)) {
    rng <- range(x, na.rm = TRUE)
    if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) == 0) {
      return(rep(mean(to), length(x)))
    }
    scales::rescale(x, to = to, from = rng)
  }

  # Normalize weights to sum to 1 if `normalize_weights` is TRUE
  weights <- c(
    directional_light_weight,
    volume_shading_weight,
    ambient_occlusion_weight
  )

  if (normalize_weights) {
    s <- sum(weights)
    if (s == 0) {
      cli::cli_abort("At least one weight must be > 0 when `normalize_weights = TRUE`.")
    }
    weights <- weights / s
  }


  # Compute directional light as the rescaled z-coordinate (light from above)
  directional_light <- safe_rescale(layout$z)

  # Compute radial volume shading as the rescaled distance from the origin
  r <- sqrt(layout$x^2 + layout$y^2 + layout$z^2)
  volume_shading <- safe_rescale(r)

  # Approximate ambient occlusion using mean distance to k nearest neighbors
  nn <- FNN::get.knn(coords, k = ambient_occlusion_k)$nn.dist
  ambient_occlusion <- safe_rescale(rowMeans(sqrt(nn)), to = c(1, 0))

  # Combine components using specified weights
  illumination <-
    weights[1] * directional_light +
    weights[2] * volume_shading +
    weights[3] * ambient_occlusion

  # Clamp illumination to specified quantiles to reduce outlier influence
  quants <- stats::quantile(illumination, probs = clamp_quantiles, na.rm = TRUE)
  illumination <- pmin(pmax(illumination, quants[[1]]), quants[[2]])

  return(illumination)
}
