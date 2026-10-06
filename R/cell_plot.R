#' Create a cell plot recipe
#'
#' Creates a cell-layout visualization recipe. The recipe stores data, column
#' mappings, and optional rendering instructions without drawing a plot. Add
#' instructions with the `cell_*()` modifier functions. Printing the recipe
#' builds and renders a ggplot. [cell_plot_interactive()] draws the same
#' recipe as a 3D Plotly widget. [cell_plot_rgl()] returns it as an rgl
#' htmlwidget. [cell_plot_animate()] encodes a rotating GIF or video after
#' [cell_coord_rotate()]. [summary()] inspects the recipe without drawing.
#'
#' A `cell_plot` is an S3 list with a fixed set of fields:
#'
#' - `data`: the visualization tibble
#' - `mapping`: named column mappings (`x`, `y`, `z`, `depth`, `color`,
#'   `size`, `alpha`, `illumination_mask`, `arrange`)
#' - `constant`: single `color`, `size`, and `alpha` values shared by every
#'   point, `NULL` when the aesthetic is mapped or left at its default
#' - `grid`, `color`, `size`, `depth`, `alpha`, `theme`, `annotation`,
#'   `illuminate`, `coord`: instruction slots, `NULL` until a modifier writes them
#'
#' Coordinate defaults are `x`, `y`, and `z`. `arrange` and `depth` default to
#' the mapped `z` column. Pass `depth = NULL` to turn off depth sizing. Color,
#' size, and alpha default to no mapping.
#'
#' `color`, `size`, and `alpha` accept either a column mapping or one constant
#' value shared by every point, such as `color = "red"`, `size = 3`, or
#' `alpha = 0.5`. A constant aesthetic has no scale, so the matching
#' `cell_node_scale_*()` modifier cannot be used with it. Without a mapping and
#' without a constant, nodes are `gray90` with relative size `1` and alpha
#' `1`.
#'
#' When a size column is mapped, the default output range is `c(2, 6)` unless
#' [cell_node_scale_size()] overrides it. For a categorical size or alpha
#' mapping, [cell_node_scale_size()] and [cell_node_scale_alpha()] also accept
#' a named vector of per-level values, such as `sizes = c(a = 1, b = 2, c = 3)`.
#' Color scale behavior is controlled by [cell_node_scale_color()]. Color scales
#' are trained once across all panels, and continuous scales span the observed
#' range unless `limits` are supplied.
#'
#' @param data A non-empty tibble containing cell-layout data.
#' @param x,y Column mappings for the horizontal and vertical coordinates.
#' Bare column names and character column names are supported.
#' @param z Column mapping for the third coordinate. Defaults to the `z`
#' column.
#' @param color Optional column mapping or constant color for every node.
#' Character, factor, numeric, and integer columns are supported, where
#' character and factor columns are treated as categorical and numeric and
#' integer columns as continuous. Missing numeric values are kept. Non-finite
#' values such as `Inf` are not. A character value is read as a column name
#' when `data` has that column and otherwise as one fully opaque color.
#' @param size,alpha Optional column mappings or constant values for node size
#' and alpha. Character, factor, numeric, and integer columns are supported,
#' where character and factor columns are treated as categorical levels.
#' Missing numeric values are kept. Non-finite values such as `Inf` are not.
#' One finite number sets a constant relative size or a constant alpha between zero
#' and one. Category-specific sizes and alphas are set with
#' [cell_node_scale_size()] and [cell_node_scale_alpha()].
#' @param arrange An optional column mapping used to order points in projected
#' plots. [cell_plot_interactive()] and [cell_plot_rgl()] ignore this mapping
#' because occlusion follows the scene camera. Defaults to the mapped `z`
#' column.
#' @param depth An optional column mapping used for ggplot depth sizing.
#' Defaults to the mapped `z` column, so depth sizing is on unless
#' `depth = NULL` turns it off. Must not be the same column as `x`, `y`, or
#' `size`. Plotly and rgl ignore this mapping.
#' @param illumination_mask Optional logical column mapping. Illumination is
#' applied only to rows where this column is `TRUE`. Rows where it is `FALSE`
#' retain their unilluminated colors. Missing values are not allowed.
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
  depth = z,
  illumination_mask = NULL
) {
  assert_class(data, "tbl_df", arg = "data")
  assert_within_limits(
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

  color_aesthetic <- .cell_plot_aesthetic(
    rlang::enquo(color),
    argument = "color",
    data = data
  )
  size_aesthetic <- .cell_plot_aesthetic(
    rlang::enquo(size),
    argument = "size",
    data = data
  )
  alpha_aesthetic <- .cell_plot_aesthetic(
    rlang::enquo(alpha),
    argument = "alpha",
    data = data
  )

  mapping <- list(
    x = .cell_plot_column(rlang::enquo(x), "x"),
    y = .cell_plot_column(rlang::enquo(y), "y"),
    z = z_mapping,
    depth = depth_mapping,
    color = color_aesthetic$column,
    size = size_aesthetic$column,
    alpha = alpha_aesthetic$column,
    illumination_mask = .cell_plot_column(
      rlang::enquo(illumination_mask),
      "illumination_mask"
    ),
    arrange = arrange_mapping
  )
  constant <- list(
    color = color_aesthetic$constant,
    size = size_aesthetic$constant,
    alpha = alpha_aesthetic$constant
  )

  mapped_columns <- unlist(mapping, use.names = FALSE)
  for (column in mapped_columns) {
    assert_col_in_data(
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
    assert_class(
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
    assert_class(
      values,
      classes = c("numeric", "integer", "character", "factor"),
      arg = column
    )
    if (is.numeric(values) && any(!is.na(values) & !is.finite(values))) {
      cli::cli_abort(
        c(
          "x" = "Mapped {.arg {aesthetic}} column {.str {column}} must contain only finite or missing values."
        )
      )
    }
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

  if (!is.null(mapping$illumination_mask)) {
    mask <- data[[mapping$illumination_mask]]
    assert_class(
      mask,
      classes = "logical",
      arg = mapping$illumination_mask
    )
    if (anyNA(mask)) {
      cli::cli_abort(
        c(
          "x" = "Mapped {.arg illumination_mask} column {.str {mapping$illumination_mask}} contains missing values."
        )
      )
    }
  }

  return(.new_cell_plot(data = data, mapping = mapping, constant = constant))
}

#' Construct a cell plot recipe
#'
#' Creates a `cell_plot` with its complete, ordered set of fields.
#'
#' @param data A non-empty tibble.
#' @param mapping A named list of normalized column mappings.
#' @param constant A named list of constant `color`, `size`, and `alpha`
#' values.
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
  constant = .cell_plot_empty_constant,
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
      constant = constant,
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

.cell_plot_empty_constant <- list(
  color = NULL,
  size = NULL,
  alpha = NULL
)

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

  assert_single_value(
    expression,
    type = "string",
    arg = argument,
    call = call
  )
  assert_non_empty_object(
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

#' Normalize a cell plot aesthetic as a column mapping or a constant
#'
#' Reads `color`, `size`, and `alpha` arguments, which accept either a column
#' name or one constant value shared by every point. Bare names are always
#' column names. A character value is a column name when `data` has that
#' column, and is otherwise read as a constant color. Numbers are constant
#' sizes and alphas.
#'
#' @param aesthetic A captured aesthetic quosure.
#' @param argument The argument name used in validation errors.
#' @param data The plot data used to tell column names from constant colors.
#' @param call The calling environment.
#'
#' @return A list with a scalar character `column` and a scalar `constant`,
#' at most one of which is non-`NULL`.
#'
#' @noRd
.cell_plot_aesthetic <- function(
  aesthetic,
  argument,
  data,
  call = rlang::caller_env()
) {
  expression <- rlang::quo_get_expr(aesthetic)
  as_column <- function() {
    list(
      column = .cell_plot_column(aesthetic, argument, call = call),
      constant = NULL
    )
  }

  if (rlang::is_null(expression) || rlang::is_symbol(expression)) {
    return(as_column())
  }
  if (rlang::is_scalar_character(expression) && expression %in% names(data)) {
    return(as_column())
  }

  value <- .cell_plot_literal(expression)
  if (identical(argument, "color")) {
    if (!rlang::is_scalar_character(value)) {
      cli::cli_abort(
        c(
          "x" = "{.arg {argument}} must be a column in {.arg data} or one color.",
          "i" = "Constant colors are single strings, as in {.code color = \"red\"}."
        ),
        call = call
      )
    }
    if (!.is_color(value)) {
      cli::cli_abort(
        c(
          "x" = "{.arg {argument}} must be a column in {.arg data} or one color.",
          "i" = "{.str {value}} is neither a column of {.arg data} nor a color."
        ),
        call = call
      )
    }
    .assert_opaque_colors(value, arg = argument, call = call)
    return(list(column = NULL, constant = value))
  }

  if (rlang::is_scalar_character(value)) {
    return(as_column())
  }
  if (!rlang::is_scalar_double(value) && !rlang::is_scalar_integer(value)) {
    cli::cli_abort(
      c(
        "x" = "{.arg {argument}} must be a column in {.arg data} or one number.",
        "i" = "Constant values are single numbers, as in {.code size = 3}."
      ),
      call = call
    )
  }

  limits <- if (identical(argument, "alpha")) c(0, 1) else c(0, Inf)
  assert_within_limits(
    value,
    limits = limits,
    arg = argument,
    call = call
  )
  if (!is.finite(value)) {
    cli::cli_abort(
      c("i" = "{.arg {argument}} must be finite."),
      call = call
    )
  }

  return(list(column = NULL, constant = as.numeric(value)))
}

#' Test whether a string names a color
#'
#' Constant colors are told apart from mistyped column names here, so an
#' unusable value can be reported as neither a column nor a color.
#'
#' @param x A scalar character value.
#'
#' @return `TRUE` when `x` can be converted to RGB.
#'
#' @noRd
.is_color <- function(x) {
  is_color <- tryCatch(
    {
      grDevices::col2rgb(x)
      TRUE
    },
    error = function(condition) FALSE
  )
  return(is_color)
}

#' Evaluate an aesthetic expression that cannot read data or variables
#'
#' Constant aesthetics are written as literals, so evaluation only needs the
#' base environment. Anything referring to a column or a variable fails here
#' and is handled as a column mapping instead.
#'
#' @param expression An aesthetic expression.
#'
#' @return The evaluated value, or `NULL` when evaluation is not possible.
#'
#' @noRd
.cell_plot_literal <- function(expression) {
  value <- tryCatch(
    eval(expression, envir = baseenv()),
    error = function(condition) NULL,
    warning = function(condition) NULL
  )
  return(value)
}

#' Validate a numeric pair used as a scale range or limits
#'
#' Checks that `x` is a length-2 finite numeric vector within `limits`, ordered
#' from low to high. This is for continuous scale output ranges (`sizes`,
#' `alphas`) and numeric scale limits, not for mapped data columns or
#' category-specific named values.
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
  assert_vector(x, type = "numeric", n = 2, arg = arg, call = call)
  assert_length(x, n = 2, arg_x = arg, call = call)
  assert_within_limits(x, limits = limits, arg = arg, call = call)
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

#' Validate a size or alpha scale output vector
#'
#' Continuous mappings require a length-2 ordered output range. Categorical
#' mappings also accept a named vector of per-level values.
#'
#' @param object A `cell_plot` recipe.
#' @param aesthetic `"size"` or `"alpha"`.
#' @param modifier The scale modifier name used in errors.
#' @param output The `sizes` or `alphas` vector.
#' @param output_arg The argument name used in validation errors.
#' @param output_limits Inclusive lower and upper bounds for `output`.
#' @param limits Optional numeric domain limits.
#' @param call The calling environment.
#'
#' @return `output`, invisibly.
#'
#' @noRd
.assert_cell_numeric_scale_output <- function(
  object,
  aesthetic,
  modifier,
  output,
  output_arg,
  output_limits,
  limits,
  call = rlang::caller_env()
) {
  .assert_cell_scale_mapping(
    object,
    aesthetic = aesthetic,
    modifier = modifier,
    call = call
  )

  values <- object$data[[object$mapping[[aesthetic]]]]
  categorical <- is.character(values) || is.factor(values)

  assert_vector(
    output,
    type = "numeric",
    n = 1,
    arg = output_arg,
    call = call
  )
  assert_within_limits(
    output,
    limits = output_limits,
    arg = output_arg,
    call = call
  )
  if (!all(is.finite(output))) {
    cli::cli_abort(
      c("i" = "{.arg {output_arg}} must contain finite values."),
      call = call
    )
  }

  if (!categorical) {
    .assert_ordered_numeric_pair(
      output,
      arg = output_arg,
      limits = output_limits,
      call = call
    )
    if (!is.null(limits)) {
      .assert_ordered_numeric_pair(
        limits,
        arg = "limits",
        limits = c(-Inf, Inf),
        call = call
      )
    }
    return(invisible(output))
  }

  if (!is.null(names(output))) {
    output_names <- names(output)
    if (any(!nzchar(output_names))) {
      cli::cli_abort(
        c(
          "i" = paste(
            "Named {.arg {output_arg}} values must name every element",
            "with a category level."
          )
        ),
        call = call
      )
    }
    assert_unique(
      output_names,
      arg = paste0("names(", output_arg, ")"),
      call = call
    )
    observed_levels <- sort(unique(as.character(values[!is.na(values)])))
    assert_x_in_y(
      observed_levels,
      output_names,
      arg_x = paste(aesthetic, "levels"),
      arg_y = paste0("names(", output_arg, ")"),
      call = call
    )
  } else {
    .assert_ordered_numeric_pair(
      output,
      arg = output_arg,
      limits = output_limits,
      call = call
    )
  }

  if (!is.null(limits)) {
    if (!is.null(names(output))) {
      cli::cli_abort(
        c(
          "x" = paste(
            "Category-specific {.arg {output_arg}} values cannot be",
            "combined with {.arg limits}."
          )
        ),
        call = call
      )
    }
    .assert_ordered_numeric_pair(
      limits,
      arg = "limits",
      limits = c(-Inf, Inf),
      call = call
    )
  }

  return(invisible(output))
}

#' Reject colors that include a translucent alpha channel
#'
#' Palette, constant, and missing-value colors must be fully opaque. Node
#' opacity belongs to `alpha` and [cell_node_scale_alpha()].
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
        "i" = "Use {.arg alpha} or {.fn cell_node_scale_alpha} to control node opacity.",
        "i" = "Translucent color(s): {.val {translucent}}"
      ),
      call = call
    )
  }

  return(invisible(x))
}

#' Require a column mapping for a cell plot scale modifier
#'
#' Scale modifiers describe how mapped values are spread across an output
#' range, so they need a mapped column. A constant aesthetic has no scale.
#'
#' @param object A `cell_plot` recipe.
#' @param aesthetic The aesthetic name.
#' @param modifier The name of the scale modifier.
#' @param call The calling environment.
#'
#' @return `object`, invisibly.
#'
#' @noRd
.assert_cell_scale_mapping <- function(
  object,
  aesthetic,
  modifier,
  call = rlang::caller_env()
) {
  if (!is.null(object$mapping[[aesthetic]])) {
    return(invisible(object))
  }

  message <- c(
    "x" = "{.fn {modifier}} requires a column mapped to {.arg {aesthetic}}."
  )
  if (!is.null(object$constant[[aesthetic]])) {
    message <- c(
      message,
      "i" = "{.arg {aesthetic}} is the constant {.val {object$constant[[aesthetic]]}}, which has no scale."
    )
  }
  cli::cli_abort(message, call = call)
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
    "data", "mapping", "constant", "grid", "color", "size", "depth", "alpha",
    "theme", "annotation", "illuminate", "coord"
  )

  assert_class(object, "cell_plot", arg = "object", call = call)
  assert_vectors_match(
    names(object),
    expected_fields,
    arg_x = "cell_plot fields",
    arg_y = "expected fields",
    call = call
  )
  assert_vectors_match(
    names(object$constant),
    names(.cell_plot_empty_constant),
    arg_x = "constant aesthetics",
    arg_y = "expected constant aesthetics",
    call = call
  )

  instruction_fields <- setdiff(
    expected_fields,
    c("data", "mapping", "constant")
  )
  for (value in object[instruction_fields]) {
    assert_class(
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

  for (aesthetic in names(.cell_plot_empty_constant)) {
    if (
      !is.null(object$constant[[aesthetic]]) &&
        !is.null(object$mapping[[aesthetic]])
    ) {
      cli::cli_abort(
        c(
          "x" = "{.arg {aesthetic}} cannot be a constant and a column mapping.",
          "i" = "Constant value: {.val {object$constant[[aesthetic]]}}",
          "i" = "Mapped column: {.str {object$mapping[[aesthetic]]}}"
        ),
        call = call
      )
    }
  }

  scale_modifiers <- c(
    color = "cell_node_scale_color",
    size = "cell_node_scale_size",
    alpha = "cell_node_scale_alpha"
  )
  for (aesthetic in names(scale_modifiers)) {
    if (!is.null(object[[aesthetic]])) {
      .assert_cell_scale_mapping(
        object,
        aesthetic = aesthetic,
        modifier = scale_modifiers[[aesthetic]],
        call = call
      )
    }
  }

  if (!is.null(object$depth) && is.null(object$mapping$depth)) {
    cli::cli_abort(
      c(
        "x" = "{.fn cell_node_depth} requires a column mapped to {.arg depth}.",
        "i" = "Depth sizing is on by default and turned off by {.code depth = NULL}."
      ),
      call = call
    )
  }
  if (!is.null(object$depth) && !is.null(object$size)) {
    cli::cli_abort(
      c(
        "x" = "{.fn cell_node_depth} and {.fn cell_node_scale_size} cannot both be used.",
        "i" = "Depth sizing scales one node size, so it has no mapped sizes to scale.",
        "i" = "Set a constant size with {.code cell_plot(size = <number>)} instead."
      ),
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
  assert_is_one_of(
    field,
    choices = setdiff(names(object), c("data", "mapping", "constant")),
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
  constants <- object$constant[
    !vapply(object$constant, is.null, logical(1))
  ]
  mapped_text <- paste(mapped, collapse = ", ")
  constants_text <- if (length(constants) == 0) {
    "none"
  } else {
    paste(names(constants), unlist(constants), sep = " = ", collapse = ", ")
  }
  modifiers_text <- if (length(modifiers) == 0) {
    "none"
  } else {
    paste(modifiers, collapse = ", ")
  }

  output <- cli::ansi_strip(c(
    cli::format_inline("{.cls cell_plot}"),
    cli::format_inline("Rows: {.val {nrow(object$data)}}"),
    cli::format_inline("Mappings: {.field {mapped_text}}"),
    cli::format_inline("Constants: {.field {constants_text}}"),
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

#' Create a complete ggplot from a `cell_plot` recipe
#'
#' @param object A `cell_plot` recipe.
#' @param ... Unused or additional arguments passed to underlying methods.
#'
#' @return A `ggplot` object.
#'
#' @export
autoplot.cell_plot <- function(object, ...) {
  plot <- object |>
    build_cell_plot() |>
    .render_cell_plot_ggplot()

  return(plot)
}
