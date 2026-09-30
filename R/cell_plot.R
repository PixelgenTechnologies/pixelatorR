#' Create a cell plot recipe
#'
#' `r lifecycle::badge("experimental")`
#'
#' Creates a cell-layout visualization recipe. The recipe stores data, column
#' mappings, and optional rendering instructions without drawing a plot. Add
#' instructions with the `cell_*()` modifier functions. Printing the recipe
#' builds and renders a ggplot. [cell_plot_interactive()] draws the same
#' recipe as a 3D Plotly widget. [cell_plot_rgl()] draws it as a native rgl
#' scene. [cell_plot_animate()] encodes a rotating GIF or video after
#' [cell_coord_rotate()]. [summary()] inspects the recipe without drawing.
#'
#' A `cell_plot` is an S3 list with a fixed set of fields:
#'
#' - `data`: the visualization tibble
#' - `mapping`: named column mappings (`x`, `y`, `z`, `depth`, `color`,
#'   `size`, `alpha`, `arrange`)
#' - `grid`, `color`, `size`, `depth`, `alpha`, `theme`, `annotation`,
#'   `illuminate`, `coord`: instruction slots, `NULL` until a modifier writes them
#'
#' Coordinate defaults are `x`, `y`, and `z`. `arrange` and `depth` default to
#' the mapped `z` column. Pass `depth = NULL` to disable depth sizing. Color,
#' size, and alpha default to no mapping.
#'
#' When no size column is mapped, nodes use a constant relative size of `1`
#' unless [cell_node_scale_size()] supplies a single `sizes` value. When a size
#' column is mapped, the default output range is `c(2, 6)` unless
#' [cell_node_scale_size()] overrides it. Color scale behavior is controlled by
#' [cell_node_scale_color()]. Color scales are trained once across all panels,
#' and continuous scales span the observed range unless `limits` are supplied.
#'
#' @param data A non-empty tibble containing cell-layout data.
#' @param x,y Column mappings for the horizontal and vertical coordinates.
#' Bare column names and character column names are supported.
#' @param z Column mapping for the third coordinate. Defaults to the `z`
#' column.
#' @param color,size,alpha Optional column mappings for node color, size, and
#' alpha. Character, factor, numeric, and integer columns are supported.
#' Character and factor columns are treated as categorical for color and as
#' discrete levels for size and alpha. Numeric and integer columns are treated
#' as continuous.
#' @param arrange An optional column mapping used to order points in projected
#' plots. [cell_plot_interactive()] and [cell_plot_rgl()] ignore this mapping
#' because occlusion follows the scene camera. Defaults to the mapped `z`
#' column.
#' @param depth An optional column mapping used for ggplot depth sizing.
#' Defaults to the mapped `z` column. Must not be the same column as `x`, `y`,
#' or `size`. Plotly and rgl ignore this mapping.
#'
#' @return A `cell_plot` recipe.
#'
#' @seealso [cell_node_scale_size()], [cell_node_depth()],
#' [cell_coord_rotate()], [cell_plot_animate()], [cell_plot_interactive()],
#' [cell_plot_rgl()]
#'
#' @examples
#' # Plot a spectral layout of a cell from the example data
#' se <- ReadPNA_Seurat(minimal_pna_pxl_file())
#' se <- LoadCellGraphs(se, cells = colnames(se)[4], verbose = FALSE) |>
#'   ComputeLayout(layout_method = "spectral")
#'
#' cell_graph <- CellGraphs(se)[[4]]
#'
#' layout_data <- FetchLayoutData(cell_graph, vars = "CD82", layout_method = "spectral_3d")
#'
#' cell_plot(layout_data, color = CD82) |>
#'   cell_node_scale_color(colors = c("lightgrey", "red")) |>
#'   cell_annotation(title = "Spectral layout", subtitle = colnames(se)[1])
#'
#' @export
cell_plot <- function(
  data,
  x = x,
  y = y,
  z = z,
  color = NULL,
  size = NULL,
  alpha = NULL,
  arrange = z,
  depth = z
) {
  pixelatorR:::assert_class(data, "tbl_df", arg = "data")
  pixelatorR:::assert_within_limits(
    nrow(data),
    limits = c(1, Inf),
    arg = "nrow(data)"
  )

  arrange_missing <- missing(arrange)
  depth_missing <- missing(depth)
  z_mapping <- .cell_plot_column(rlang::enquo(z), "z")

  arrange_mapping <- if (arrange_missing) {
    z_mapping
  } else {
    .cell_plot_column(rlang::enquo(arrange), "arrange")
  }

  depth_mapping <- if (depth_missing) {
    z_mapping
  } else {
    .cell_plot_column(rlang::enquo(depth), "depth")
  }

  mapping <- list(
    x = .cell_plot_column(rlang::enquo(x), "x"),
    y = .cell_plot_column(rlang::enquo(y), "y"),
    z = z_mapping,
    depth = depth_mapping,
    color = .cell_plot_column(rlang::enquo(color), "color"),
    size = .cell_plot_column(rlang::enquo(size), "size"),
    alpha = .cell_plot_column(rlang::enquo(alpha), "alpha"),
    arrange = arrange_mapping
  )

  mapped_columns <- unlist(mapping, use.names = FALSE)
  for (column in mapped_columns) {
    pixelatorR:::assert_col_in_data(
      column,
      data = data,
      arg_data = "data"
    )
  }

  coordinate_columns <- unlist(
    mapping[c("x", "y", "z", "depth")],
    use.names = FALSE
  )
  for (column in coordinate_columns) {
    pixelatorR:::assert_class(
      data[[column]],
      classes = c("numeric", "integer"),
      arg = column
    )
    if (!all(is.finite(data[[column]]))) {
      cli::cli_abort(
        c(
          "x" = "Coordinate column {.str {column}} must contain only finite values."
        )
      )
    }
  }

  for (aesthetic in c("color", "size", "alpha")) {
    column <- mapping[[aesthetic]]
    if (is.null(column)) {
      next
    }
    values <- data[[column]]
    pixelatorR:::assert_class(
      values,
      classes = c("numeric", "integer", "character", "factor"),
      arg = column
    )
    has_usable_value <- if (is.numeric(values)) {
      any(is.finite(values))
    } else {
      any(!is.na(values))
    }
    if (!has_usable_value) {
      cli::cli_abort(
        c("x" = "Mapped {.arg {aesthetic}} column {.str {column}} has no usable values.")
      )
    }
  }

  return(.new_cell_plot(data = data, mapping = mapping))
}

#' Construct a cell plot recipe
#'
#' Creates a `cell_plot` with its complete, ordered set of fields.
#'
#' @param data A non-empty tibble.
#' @param mapping A named list of normalized column mappings.
#' @param grid,color,size,depth,alpha,theme,annotation,illuminate,coord Optional
#' instruction specifications.
#' @param call The calling environment.
#'
#' @return A validated `cell_plot` object.
#'
#' @noRd
.new_cell_plot <- function(
  data,
  mapping,
  grid = NULL,
  color = NULL,
  size = NULL,
  depth = NULL,
  alpha = NULL,
  theme = NULL,
  annotation = NULL,
  illuminate = NULL,
  coord = NULL,
  call = rlang::caller_env()
) {
  object <- structure(
    list(
      data = data,
      mapping = mapping,
      grid = grid,
      color = color,
      size = size,
      depth = depth,
      alpha = alpha,
      theme = theme,
      annotation = annotation,
      illuminate = illuminate,
      coord = coord
    ),
    class = "cell_plot"
  )

  .validate_cell_plot(object, call = call)
  return(object)
}

.cell_plot_default_theme <- list(
  background_color = "white",
  strip_background_color = "#D9D9D9",
  text_color = "black",
  text_size = 11
)

# Match the explicit spacing and row-label direction used by ggplot2's
# default facet strips across cell plot renderers.
.cell_plot_facet_strip_margin <- ggplot2::margin(4.4, 4.4, 4.4, 4.4)
.cell_plot_row_strip_angle <- 90

.cell_plot_default_annotation <- list(
  title = NULL,
  subtitle = NULL,
  legend_title = NULL
)

#' Resolve the color legend title for a cell plot
#'
#' Uses [cell_annotation()] `legend_title` when it is set. Otherwise the mapped
#' color column name is used.
#'
#' @param object A `cell_plot` or `cell_plot_built` object.
#'
#' @return A scalar character legend title, or `NULL` when no color mapping
#' exists.
#'
#' @noRd
.cell_plot_legend_title <- function(object) {
  object$annotation$legend_title %||% object$mapping$color
}

#' Normalize a cell plot column mapping
#'
#' Converts a bare or character column name to a scalar character value.
#'
#' @param mapping A captured mapping quosure.
#' @param argument The argument name used in validation errors.
#' @param call The calling environment.
#'
#' @return A scalar character column name or `NULL`.
#'
#' @noRd
.cell_plot_column <- function(mapping, argument, call = rlang::caller_env()) {
  expression <- rlang::quo_get_expr(mapping)

  if (rlang::is_null(expression)) {
    return(NULL)
  }
  if (rlang::is_symbol(expression)) {
    return(rlang::as_string(expression))
  }

  pixelatorR:::assert_single_value(
    expression,
    type = "string",
    arg = argument,
    call = call
  )
  pixelatorR:::assert_non_empty_object(
    expression,
    classes = "character",
    arg = argument,
    call = call
  )

  if (!nzchar(expression)) {
    cli::cli_abort(
      c("i" = "{.arg {argument}} must not be an empty column name."),
      call = call
    )
  }

  return(expression)
}

#' Validate a numeric pair used as a scale range or limits
#'
#' Checks that `x` is a length-2 finite numeric vector within `limits`, ordered
#' from low to high. This is for scale output ranges (`sizes`, `alphas`) and
#' numeric scale limits, not for mapped data columns.
#'
#' @param x A numeric vector of length 2.
#' @param arg The argument name used in validation errors.
#' @param limits Inclusive lower and upper bounds.
#' @param strictly_increasing Whether the first value must be smaller than the second.
#' @param call The calling environment.
#'
#' @return `x`, invisibly.
#'
#' @noRd
.assert_ordered_numeric_pair <- function(
  x,
  arg,
  limits,
  strictly_increasing = FALSE,
  call = rlang::caller_env()
) {
  pixelatorR:::assert_vector(x, type = "numeric", n = 2, arg = arg, call = call)
  pixelatorR:::assert_length(x, n = 2, arg_x = arg, call = call)
  pixelatorR:::assert_within_limits(x, limits = limits, arg = arg, call = call)
  if (!all(is.finite(x))) {
    cli::cli_abort(
      c("i" = "{.arg {arg}} must contain finite values."),
      call = call
    )
  }
  if (strictly_increasing && x[1] >= x[2]) {
    cli::cli_abort(
      c("i" = "{.arg {arg}} must contain two increasing values."),
      call = call
    )
  }
  if (!strictly_increasing && x[1] > x[2]) {
    cli::cli_abort(
      c("i" = "{.arg {arg}} must be ordered from low to high."),
      call = call
    )
  }

  return(invisible(x))
}

#' Reject colors that include a translucent alpha channel
#'
#' Palette and missing-value colors must be fully opaque. Node opacity belongs
#' to the alpha mapping and [cell_node_scale_alpha()].
#'
#' @param x A vector of valid colors.
#' @param arg The argument name used in validation errors.
#' @param call The calling environment.
#'
#' @return `x`, invisibly.
#'
#' @noRd
.assert_opaque_colors <- function(x, arg, call = rlang::caller_env()) {
  alpha <- grDevices::col2rgb(x, alpha = TRUE)[4, ]
  translucent <- x[alpha != 255]
  if (length(translucent) > 0) {
    cli::cli_abort(
      c(
        "x" = "{.arg {arg}} must be fully opaque.",
        "i" = "Use {.fn cell_node_scale_alpha} to control node opacity.",
        "i" = "Translucent color(s): {.val {translucent}}"
      ),
      call = call
    )
  }

  return(invisible(x))
}

#' Validate a cell plot recipe
#'
#' Checks class, field names, instruction types, and mapping dependencies.
#' Data and mapping contents are validated when they are first set.
#'
#' @param object An object to validate.
#' @param call The calling environment.
#'
#' @return `object`, invisibly.
#'
#' @noRd
.validate_cell_plot <- function(object, call = rlang::caller_env()) {
  expected_fields <- c(
    "data", "mapping", "grid", "color", "size", "depth", "alpha", "theme",
    "annotation", "illuminate", "coord"
  )

  pixelatorR:::assert_class(object, "cell_plot", arg = "object", call = call)
  pixelatorR:::assert_vectors_match(
    names(object),
    expected_fields,
    arg_x = "cell_plot fields",
    arg_y = "expected fields",
    call = call
  )

  instruction_fields <- setdiff(expected_fields, c("data", "mapping"))
  for (value in object[instruction_fields]) {
    pixelatorR:::assert_class(
      value,
      classes = "list",
      allow_null = TRUE,
      arg = "instruction",
      call = call
    )
  }

  if (
    is.null(object$mapping$x) ||
      is.null(object$mapping$y) ||
      is.null(object$mapping$z)
  ) {
    cli::cli_abort(
      c("x" = "{.arg x}, {.arg y}, and {.arg z} mappings are required."),
      call = call
    )
  }

  if (!is.null(object$color) && is.null(object$mapping$color)) {
    cli::cli_abort(
      c("x" = "{.fn cell_node_scale_color} requires a column mapped to {.arg color}."),
      call = call
    )
  }
  if (!is.null(object$size)) {
    sizes <- object$size$sizes
    mapped_size <- !is.null(object$mapping$size)
    if (length(sizes) == 1L) {
      if (mapped_size) {
        cli::cli_abort(
          c(
            "x" = paste0(
              "{.fn cell_node_scale_size} with a single {.arg sizes} value ",
              "cannot be used when a column is mapped to {.arg size}."
            )
          ),
          call = call
        )
      }
      if (!is.null(object$size$limits)) {
        cli::cli_abort(
          c("x" = "{.arg limits} is only used with a size mapping."),
          call = call
        )
      }
    } else if (!mapped_size) {
      cli::cli_abort(
        c("x" = "{.fn cell_node_scale_size} requires a column mapped to {.arg size}."),
        call = call
      )
    }
  }
  if (!is.null(object$depth) && is.null(object$mapping$depth)) {
    cli::cli_abort(
      c("x" = "{.fn cell_node_depth} requires a column mapped to {.arg depth}."),
      call = call
    )
  }
  if (!is.null(object$mapping$depth)) {
    for (aesthetic in c("x", "y", "size")) {
      mapped_column <- object$mapping[[aesthetic]]
      if (
        !is.null(mapped_column) &&
          identical(object$mapping$depth, mapped_column)
      ) {
        cli::cli_abort(
          c(
            "x" = "{.arg depth} cannot use the same column as {.arg {aesthetic}}.",
            "i" = "Both are mapped to {.str {object$mapping$depth}}."
          ),
          call = call
        )
      }
    }
  }
  if (!is.null(object$alpha) && is.null(object$mapping$alpha)) {
    cli::cli_abort(
      c("x" = "{.fn cell_node_scale_alpha} requires a column mapped to {.arg alpha}."),
      call = call
    )
  }
  if (!is.null(object$coord) && is.null(object$mapping$z)) {
    cli::cli_abort(
      c("x" = "{.fn cell_coord_rotate} requires a {.arg z} coordinate."),
      call = call
    )
  }

  return(invisible(object))
}

#' Replace one cell plot instruction
#'
#' Returns a new recipe after replacing one existing instruction field.
#'
#' @param object A `cell_plot` recipe.
#' @param field The instruction field to replace.
#' @param value A list specification or `NULL`.
#' @param call The calling environment.
#'
#' @return A new `cell_plot` recipe.
#'
#' @noRd
.replace_cell_plot_field <- function(
  object,
  field,
  value,
  call = rlang::caller_env()
) {
  pixelatorR:::assert_is_one_of(
    field,
    choices = setdiff(names(object), c("data", "mapping")),
    arg = "field",
    call = call
  )

  object[[field]] <- value
  .validate_cell_plot(object, call = call)
  return(object)
}

#' Summarise a cell plot recipe
#'
#' Creates a compact summary of a recipe without building or drawing it.
#'
#' @param object A `cell_plot` recipe.
#' @param ... Additional arguments. Currently not used.
#'
#' @return A `summary.cell_plot` object. Printing writes one line per
#' field. The underlying value is a character vector describing the recipe.
#'
#' @export
summary.cell_plot <- function(object, ...) {
  .validate_cell_plot(object)

  mapped <- names(object$mapping)[!vapply(object$mapping, is.null, logical(1))]
  modifiers <- names(object)[
    names(object) %in% c(
      "grid", "color", "size", "depth", "alpha", "theme", "annotation",
      "illuminate", "coord"
    ) &
      !vapply(object, is.null, logical(1))
  ]
  mapped_text <- paste(mapped, collapse = ", ")
  modifiers_text <- if (length(modifiers) == 0) {
    "none"
  } else {
    paste(modifiers, collapse = ", ")
  }

  output <- cli::ansi_strip(c(
    cli::format_inline("{.cls cell_plot}"),
    cli::format_inline("Rows: {.val {nrow(object$data)}}"),
    cli::format_inline("Mappings: {.field {mapped_text}}"),
    cli::format_inline("Modifiers: {.field {modifiers_text}}")
  ))

  return(structure(output, class = c("summary.cell_plot", "character")))
}

#' Print a cell plot summary
#'
#' Writes the recipe summary with one field per line.
#'
#' @param x A `summary.cell_plot` object.
#' @param ... Additional arguments. Currently not used.
#'
#' @return `x`, invisibly.
#'
#' @export
print.summary.cell_plot <- function(x, ...) {
  writeLines(x)
  return(invisible(x))
}

#' Print a cell plot recipe
#'
#' Builds and draws the recipe as a static ggplot, then returns that ggplot
#' invisibly. This matches `print.ggplot()`: typing a recipe at the prompt
#' renders the plot. Use [cell_plot_interactive()] for the Plotly renderer or
#' [cell_plot_rgl()] for the native rgl renderer.
#'
#' @param x A `cell_plot` recipe.
#' @param ... Additional arguments. Currently not used.
#'
#' @return A `ggplot` object, invisibly.
#'
#' @export
print.cell_plot <- function(x, ...) {
  plot <- x |>
    build_cell_plot() |>
    .render_cell_plot_ggplot()
  return(print(plot))
}
