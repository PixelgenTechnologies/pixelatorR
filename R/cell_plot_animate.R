#' Render a rotating cell plot animation
#'
#' `r lifecycle::badge("experimental")`
#'
#' Builds a [cell_plot()] recipe once, draws one frame per rotation angle, and
#' encodes a GIF or video. Rotation geometry comes from [cell_coord_rotate()].
#' File type, size, resolution, frame rate, and the frame backend belong here.
#'
#' GIF output uses gifski. Any other extension is encoded with av. The parent
#' directory of `file` must already exist. An existing file is overwritten.
#'
#' @param object A `cell_plot` recipe that includes [cell_coord_rotate()].
#' @param file Output path. The extension selects the encoder.
#' @param frames Positive whole number of encoded frames.
#' @param width,height Output size in pixels.
#' @param res PNG resolution in pixels per inch.
#' @param fps Encoded frames per second.
#' @param frame_backend Device used to draw each frame. `"base"` is faster.
#' `"ggplot2"` matches the static renderer more closely.
#' @param workers Positive whole number of parallel workers. `1` renders
#' sequentially.
#'
#' @return The output path, invisibly.
#'
#' @seealso [cell_plot()], [cell_coord_rotate()], [cell_illuminate()]
#'
#' @examplesIf requireNamespace("gifski", quietly = TRUE)
#' se <- ReadPNA_Seurat(minimal_pna_pxl_file())
#' se <- LoadCellGraphs(se, cells = colnames(se)[4], verbose = FALSE) |>
#'   ComputeLayout(layout_method = "spectral")
#'
#' cell_graph <- CellGraphs(se)[[4]]
#'
#' layout_data <- FetchLayoutData(cell_graph, layout_method = "spectral_3d") |>
#'   # Downsample to speed up tests
#'   dplyr::slice_sample(n = 5000)
#' file <- tempfile(fileext = ".gif")
#' cell_plot(layout_data) |>
#'   cell_coord_rotate(axis = "y") |>
#'   cell_plot_animate(file, frames = 2, width = 160, height = 160, res = 72)
#'
#' @export
cell_plot_animate <- function(
  object,
  file,
  frames = 500,
  width = 500,
  height = 500,
  res = 150,
  fps = 20,
  frame_backend = c("base", "ggplot2"),
  workers = 1L
) {
  .validate_cell_plot(object)
  if (is.null(object$coord) || !identical(object$coord$type, "rotate")) {
    cli::cli_abort(
      c(
        "x" = "{.fn cell_plot_animate} requires {.fn cell_coord_rotate}.",
        "i" = "Add rotation geometry before animating."
      )
    )
  }

  pixelatorR:::assert_single_value(file, type = "string", arg = "file")
  if (!nzchar(fs::path_ext(file))) {
    cli::cli_abort(
      c("x" = "{.arg file} must include an extension such as {.str gif} or {.str mp4}.")
    )
  }
  parent_dir <- fs::path_dir(fs::path_abs(file))
  if (!fs::dir_exists(parent_dir)) {
    cli::cli_abort(
      c(
        "x" = "The parent directory of {.arg file} does not exist.",
        "i" = "Missing directory: {.path {parent_dir}}"
      )
    )
  }
  for (argument in c("frames", "width", "height", "workers")) {
    pixelatorR:::assert_single_value(
      get(argument),
      type = "integer",
      arg = argument
    )
    pixelatorR:::assert_within_limits(
      get(argument),
      limits = c(1, Inf),
      arg = argument
    )
  }
  for (argument in c("res", "fps")) {
    pixelatorR:::assert_single_value(
      get(argument),
      type = "numeric",
      arg = argument
    )
    pixelatorR:::assert_within_limits(
      get(argument),
      limits = c(.Machine$double.eps, Inf),
      arg = argument
    )
    if (!is.finite(get(argument))) {
      cli::cli_abort(c("i" = "{.arg {argument}} must be finite."))
    }
  }
  frame_backend <- match.arg(frame_backend)
  frames <- as.integer(frames)
  width <- as.integer(width)
  height <- as.integer(height)
  workers <- as.integer(workers)

  extension <- tolower(fs::path_ext(file))
  if (identical(extension, "gif")) {
    rlang::check_installed("gifski")
  } else {
    rlang::check_installed("av")
  }

  built <- build_cell_plot(object)
  built <- .cell_prepare_animation_illumination(built)
  angles <- .cell_animation_angles(built$coord, frames = frames)
  limits <- .cell_animation_plot_limits(.cell_animation_limits(built, angles))
  png_device <- .cell_animation_png_device()
  tmp_dir <- fs::file_temp("cell_plot_frames")
  fs::dir_create(tmp_dir)
  png_files <- fs::path(
    tmp_dir,
    sprintf("frame_%04d.png", seq_along(angles))
  )
  cluster <- NULL
  encoded <- FALSE
  on.exit(
    {
      if (!is.null(cluster)) {
        try(parallel::stopCluster(cluster), silent = TRUE)
      }
    },
    add = TRUE
  )

  tryCatch(
    {
      cli::cli_progress_bar("Rendering frames", total = length(angles))
      if (workers == 1L) {
        for (i in seq_along(angles)) {
          .cell_animation_write_frame(
            object = built,
            angle = angles[[i]],
            path = png_files[[i]],
            limits = limits,
            width = width,
            height = height,
            res = res,
            frame_backend = frame_backend,
            png_device = png_device
          )
          cli::cli_progress_update()
        }
      } else {
        workers <- min(workers, length(angles))
        cluster <- parallel::makeCluster(workers)
        parallel::clusterEvalQ(cluster, {
          for (package in c(
            "grDevices",
            "graphics",
            "ggplot2",
            "tibble",
            "dplyr",
            "scales",
            "cli",
            "rlang",
            "fs",
            "pixelatorR"
          )) {
            library(package, character.only = TRUE)
          }
          if (requireNamespace("ragg", quietly = TRUE)) {
            library("ragg", character.only = TRUE)
          }
          NULL
        })
        worker_fun <- .cell_animation_parallel_fun(
          built = built,
          angles = angles,
          png_files = as.character(png_files),
          limits = limits,
          width = width,
          height = height,
          res = res,
          frame_backend = frame_backend,
          png_device = png_device
        )
        parallel::parLapplyLB(cluster, seq_along(angles), worker_fun)
        cli::cli_progress_update(inc = length(angles))
        parallel::stopCluster(cluster)
        cluster <- NULL
      }
      cli::cli_progress_done()

      cli::cli_progress_step("Encoding {.file {file}}")
      if (identical(extension, "gif")) {
        gifski::gifski(
          png_files = png_files,
          gif_file = file,
          width = width,
          height = height,
          delay = 1 / fps
        )
      } else {
        av::av_encode_video(
          input = png_files,
          output = file,
          framerate = fps
        )
      }
      encoded <- TRUE
    },
    error = function(cnd) {
      cli::cli_abort(
        c(
          "x" = "Failed to render or encode the animation.",
          "i" = "Rendered frames were kept in {.path {tmp_dir}}."
        ),
        parent = cnd
      )
    }
  )

  if (encoded) {
    try(fs::dir_delete(tmp_dir), silent = TRUE)
  }
  return(invisible(file))
}

#' Choose a PNG device for animation frames
#'
#' Prefers ragg when it is installed.
#'
#' @return A PNG device function.
#'
#' @noRd
.cell_animation_png_device <- function() {
  if (requireNamespace("ragg", quietly = TRUE)) {
    return(ragg::agg_png)
  }
  cli::cli_inform(
    "For faster rendering, install the {.pkg ragg} package."
  )
  return(grDevices::png)
}

#' Expand degenerate projected animation limits
#'
#' Keeps a shared viewport finite when a rotated axis has no range.
#'
#' @param limits A list with numeric `x` and `y` ranges.
#'
#' @return The same list with finite, non-zero ranges.
#'
#' @noRd
.cell_animation_plot_limits <- function(limits) {
  pad_range <- function(rng) {
    if (length(rng) != 2L || !all(is.finite(rng))) {
      return(c(-1, 1))
    }
    if (diff(rng) == 0) {
      pad <- max(abs(rng[1]), 1) * 0.05
      return(rng + c(-pad, pad))
    }
    pad <- diff(rng) * 0.05
    return(rng + c(-pad, pad))
  }
  return(list(x = pad_range(limits$x), y = pad_range(limits$y)))
}

#' Write one animation frame to a PNG file
#'
#' @param object A `cell_plot_built` object.
#' @param angle Frame angle in degrees.
#' @param path Output PNG path.
#' @param limits Shared plot limits.
#' @param width,height,res PNG device settings.
#' @param frame_backend `"base"` or `"ggplot2"`.
#' @param png_device PNG device function.
#'
#' @return `path`, invisibly.
#'
#' @noRd
.cell_animation_write_frame <- function(
  object,
  angle,
  path,
  limits,
  width,
  height,
  res,
  frame_backend,
  png_device
) {
  frame <- .cell_animation_frame(object, angle)
  background <- frame$theme$background_color
  png_device(
    filename = path,
    width = width,
    height = height,
    res = res,
    units = "px",
    bg = background
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  if (identical(frame_backend, "ggplot2")) {
    print(.render_cell_plot_ggplot(frame, limits = limits))
  } else {
    .render_cell_plot_base(frame, limits = limits)
  }
  return(invisible(path))
}

#' Build a self-contained worker that writes one animation frame
#'
#' Copies animation helpers into a plain environment so PSOCK workers do not
#' depend on the package namespace.
#'
#' @param built A `cell_plot_built` object.
#' @param angles Frame angles in degrees.
#' @param png_files Output PNG paths.
#' @param limits Shared plot limits.
#' @param width,height,res PNG device settings.
#' @param frame_backend `"base"` or `"ggplot2"`.
#' @param png_device PNG device function.
#'
#' @return A function of a frame index.
#'
#' @noRd
.cell_animation_parallel_fun <- function(
  built,
  angles,
  png_files,
  limits,
  width,
  height,
  res,
  frame_backend,
  png_device
) {
  namespace <- asNamespace("pixelatorR")
  export_names <- grep(
    "^\\.(cell_|normalize_cell_direction|render_cell_plot_)",
    ls(namespace, all.names = TRUE),
    value = TRUE
  )
  worker_env <- new.env(parent = globalenv())
  worker_env$built <- built
  worker_env$angles <- angles
  worker_env$png_files <- png_files
  worker_env$limits <- limits
  worker_env$width <- width
  worker_env$height <- height
  worker_env$res <- res
  worker_env$frame_backend <- frame_backend
  worker_env$png_device <- png_device
  list2env(mget(export_names, envir = namespace), envir = worker_env)
  for (name in export_names) {
    if (is.function(worker_env[[name]])) {
      environment(worker_env[[name]]) <- worker_env
    }
  }
  fun <- function(i) {
    .cell_animation_write_frame(
      object = built,
      angle = angles[[i]],
      path = png_files[[i]],
      limits = limits,
      width = width,
      height = height,
      res = res,
      frame_backend = frame_backend,
      png_device = png_device
    )
  }
  environment(fun) <- worker_env
  return(fun)
}

#' Create frame angles for a cell plot rotation
#'
#' Generates a forward or boomerang sequence with exactly `frames` entries.
#' Full rotations omit the duplicate closing angle.
#'
#' @param specification A rotation specification from [cell_coord_rotate()].
#' @param frames Positive whole number of output frames.
#'
#' @return A numeric vector of angles in degrees.
#'
#' @noRd
.cell_animation_angles <- function(specification, frames) {
  pixelatorR:::assert_single_value(
    frames,
    type = "integer",
    arg = "frames"
  )
  pixelatorR:::assert_within_limits(
    frames,
    limits = c(1, Inf),
    arg = "frames"
  )
  frames <- as.integer(frames)
  max_degree <- specification$max_degree

  if (frames == 1L) {
    return(0)
  }

  if (!specification$boomerang) {
    angles <- if (abs(max_degree) == 360) {
      seq(0, max_degree, length.out = frames + 1L)[-(frames + 1L)]
    } else {
      seq(0, max_degree, length.out = frames)
    }
    return(angles)
  }

  outward_frames <- floor(frames / 2) + 1L
  angles <- seq(0, max_degree, length.out = outward_frames)
  return_frames <- frames - outward_frames
  if (return_frames == 0L) {
    return(angles)
  }

  return_angles <- seq(
    max_degree,
    0,
    length.out = return_frames + 2L
  )
  return_angles <- return_angles[-c(1L, length(return_angles))]
  return(c(angles, return_angles))
}

#' Convert a cell plot rotation axis to a unit vector
#'
#' Expands a principal-axis name or normalizes an arbitrary numeric axis.
#'
#' @param axis `"x"`, `"y"`, `"z"`, or a numeric direction of length three.
#'
#' @return A numeric unit vector in x, y, and z order.
#'
#' @noRd
.cell_rotation_axis <- function(axis) {
  if (is.character(axis)) {
    direction <- switch(axis,
      x = c(1, 0, 0),
      y = c(0, 1, 0),
      z = c(0, 0, 1)
    )
    return(direction)
  }

  return(.normalize_cell_direction(axis, argument = "axis"))
}

#' Rotate three-dimensional cell plot coordinates
#'
#' Applies an axis-angle rotation around the coordinate origin or independently
#' around each panel centroid.
#'
#' @param layout A tibble with numeric `x`, `y`, and `z` columns.
#' @param axis Principal-axis name or numeric direction.
#' @param angle Angle in degrees following the right-hand rule.
#' @param origin `"origo"` or `"centroid"`.
#' @param panel_rows A list of integer row indices for non-empty panels.
#'
#' @return A tibble of rotated `x`, `y`, and `z` coordinates.
#'
#' @noRd
.cell_rotate_layout <- function(
  layout,
  axis,
  angle,
  origin,
  panel_rows
) {
  axis <- .cell_rotation_axis(axis)
  angle <- angle * pi / 180
  cross_product <- matrix(
    c(
      0, -axis[3], axis[2],
      axis[3], 0, -axis[1],
      -axis[2], axis[1], 0
    ),
    nrow = 3,
    byrow = TRUE
  )
  rotation <- diag(3) * cos(angle) +
    (1 - cos(angle)) * tcrossprod(axis) +
    sin(angle) * cross_product

  coordinates <- as.matrix(layout[, c("x", "y", "z")])
  for (rows in panel_rows) {
    pivot <- if (origin == "centroid") {
      colMeans(coordinates[rows, , drop = FALSE])
    } else {
      c(0, 0, 0)
    }
    centered <- sweep(
      coordinates[rows, , drop = FALSE],
      MARGIN = 2,
      STATS = pivot
    )
    coordinates[rows, ] <- sweep(
      centered %*% t(rotation),
      MARGIN = 2,
      STATS = pivot,
      FUN = "+"
    )
  }

  return(tibble::as_tibble(
    coordinates,
    .name_repair = ~ c("x", "y", "z")
  ))
}

#' Find fixed projected limits for a cell plot animation
#'
#' Rotates coordinates through the requested frame sequence and finds one
#' shared x and y range without changing the coordinate scale.
#'
#' @param object A `cell_plot_built` object with a rotation specification.
#' @param angles Numeric frame angles in degrees.
#'
#' @return A list containing numeric `x` and `y` ranges.
#'
#' @noRd
.cell_animation_limits <- function(object, angles) {
  pixelatorR:::assert_class(object, "cell_plot_built", arg = "object")
  mapping <- object$mapping
  layout <- tibble::tibble(
    x = object$data[[mapping$x]],
    y = object$data[[mapping$y]],
    z = object$data[[mapping$z]]
  )
  panel_rows <- .cell_plot_panel_rows(object$data, object$grid)
  limits <- list(x = numeric(), y = numeric())

  for (angle in angles) {
    rotated <- .cell_rotate_layout(
      layout = layout,
      axis = object$coord$axis,
      angle = angle,
      origin = object$coord$origin,
      panel_rows = panel_rows
    )
    limits$x <- range(limits$x, rotated$x)
    limits$y <- range(limits$y, rotated$y)
  }

  return(limits)
}

#' Reorder row-aligned values in a built scale
#'
#' Reorders scale fields that contain one value per data row while preserving
#' scalar defaults and scale metadata.
#'
#' @param scale A built color, size, or alpha scale.
#' @param order Integer row order.
#' @param rows Number of rows in the built plot.
#'
#' @return The reordered scale.
#'
#' @noRd
.cell_reorder_scale <- function(scale, order, rows) {
  for (field in c("values", "resolved", "illuminated")) {
    if (!is.null(scale[[field]]) && length(scale[[field]]) == rows) {
      scale[[field]] <- scale[[field]][order]
    }
  }
  return(scale)
}

#' Create collision-free columns for rotated coordinates
#'
#' Generates internal column names that do not overwrite columns supplied by
#' the user. These columns keep rotated x, y, and z values independent when
#' multiple coordinate mappings refer to the same source column.
#'
#' @param data A data frame containing the plot data.
#'
#' @return Three character column names for x, y, and z.
#'
#' @noRd
.cell_animation_coordinate_names <- function(data) {
  coordinate_names <- paste0(".cell_animation_", c("x", "y", "z"))
  while (any(coordinate_names %in% names(data))) {
    coordinate_names <- paste0(coordinate_names, "_")
  }
  return(coordinate_names)
}

#' Prepare one rotated cell plot animation frame
#'
#' Rotates coordinates, optionally recomputes locked-light illumination, sorts
#' points by camera depth, and keeps row-aligned aesthetics synchronized.
#'
#' @param object A `cell_plot_built` object with a rotation specification.
#' @param angle Frame angle in degrees.
#'
#' @return A `cell_plot_built` object prepared for frame rendering.
#'
#' @noRd
.cell_animation_frame <- function(object, angle) {
  pixelatorR:::assert_class(object, "cell_plot_built", arg = "object")
  object <- .cell_prepare_animation_illumination(object)
  illumination_cache <- attr(
    object,
    "cell_illumination_cache",
    exact = TRUE
  )
  mapping <- object$mapping
  rows <- nrow(object$data)
  panel_rows <- .cell_plot_panel_rows(object$data, object$grid)
  layout <- tibble::tibble(
    x = object$data[[mapping$x]],
    y = object$data[[mapping$y]],
    z = object$data[[mapping$z]]
  )
  rotated <- .cell_rotate_layout(
    layout = layout,
    axis = object$coord$axis,
    angle = angle,
    origin = object$coord$origin,
    panel_rows = panel_rows
  )

  coordinate_mappings <- unname(unlist(mapping[c("x", "y", "z")]))
  if (anyDuplicated(coordinate_mappings)) {
    coordinate_mappings <- .cell_animation_coordinate_names(object$data)
    object$mapping$x <- coordinate_mappings[[1]]
    object$mapping$y <- coordinate_mappings[[2]]
    object$mapping$z <- coordinate_mappings[[3]]
  }
  object$data[[coordinate_mappings[[1]]]] <- rotated$x
  object$data[[coordinate_mappings[[2]]]] <- rotated$y
  object$data[[coordinate_mappings[[3]]]] <- rotated$z

  if (!is.null(object$illuminate) && isTRUE(object$illuminate$lock_light)) {
    object$color$illuminated <- .cell_apply_cached_illumination(
      layout = rotated,
      colors = object$color$resolved,
      specification = object$illuminate,
      cache = illumination_cache
    )
  }

  depth_order <- order(rotated$z, na.last = TRUE)
  object$data <- object$data[depth_order, , drop = FALSE]
  object$color <- .cell_reorder_scale(object$color, depth_order, rows)
  object$size <- .cell_reorder_scale(object$size, depth_order, rows)
  object$alpha <- .cell_reorder_scale(object$alpha, depth_order, rows)

  if (!is.null(mapping$depth)) {
    object$mapping$depth <- object$mapping$z
  }

  attr(object, "cell_illumination_cache") <- NULL
  return(object)
}
