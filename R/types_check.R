#' Type check helpers
#'
#' Utility functions for type checking in \code{pixelatorR}. These functions
#' throw informative error messages if the checks fail.
#'
#' @param x,y An object to check
#' @param data A data-frame like object
#' @param type A character string specifying the type or vector class to check against.
#' @param allow_null Either \code{TRUE} or \code{FALSE}. If \code{TRUE}, \code{x}
#' can be \code{NULL} and the check passes.
#' @param n An integer
#' @param arg,arg_x,arg_y,arg_data The name of an argument to check. Used for error messages.
#' @param call An environment, typically the environment in which the function
#' was called
#' @param classes A character vector of classes to check against
#' @param limits A numeric vector of length 2 specifying a lower and upper limit
#' @param ext A character string specifying a file extension
#'
#' @returns Nothing if the check passes, otherwise throws an error.
#'
#' @details
#' - \code{assert_single_value} checks if \code{x} is a single value of a specified \code{type}.
#' - \code{assert_vector} checks if \code{x} is a vector of a specified \code{type} and with at
#' least \code{n} elements.
#' - \code{assert_class} checks if \code{x} is of a specific class.
#' - \code{assert_mpx_assay} checks if \code{x} is a \code{CellGraphAssay} or \code{CellGraphAssay5}.
#' - \code{assert_x_in_y} checks if all elements of \code{x} are in \code{y}.
#' - \code{assert_single_values_are_different} checks if \code{x} and \code{y} are not different strings.
#' - \code{assert_col_class} checks if column \code{x} in \code{data} is of a specific class.
#' - \code{assert_col_in_data} checks if column \code{x} is present in \code{data}.
#' - \code{assert_vectors_match} checks if vectors \code{x} and \code{y} are identical.
#' - \code{assert_length} checks if vector \code{x} has a specific length.
#' - \code{assert_within_limits} checks if values in \code{x} are within \code{limits}.
#' - \code{assert_function} checks \code{x} is a function.
#' - \code{assert_file_exists} checks is file \code{x} exists.
#' - \code{assert_file_ext} checks is file \code{x} has file extension \code{ext}.
#' - \code{assert_unique} checks if values in \code{x} are unique.
#'
#' @name type-checks
#' @rdname type_check_helpers
#'
assert_single_value <- function(
  x,
  type = c("string", "numeric", "integer", "bool"),
  allow_null = FALSE,
  arg = caller_arg(x),
  call = caller_env()
) {
  if (allow_null && is.null(x)) {
    return(invisible(NULL))
  }
  type <- match.arg(type, choices = c("string", "numeric", "integer", "bool"))
  check <- switch(type,
    "string" = \(x) {
      rlang::is_string(x)
    },
    "numeric" = \(x) {
      is.numeric(x) && length(x) == 1
    },
    "integer" = \(x) {
      is_scalar_wholenumber(x)
    },
    "bool" = \(x) {
      rlang::is_scalar_logical(x)
    }
  )
  msg <- switch(type,
    "string" = "a single string",
    "numeric" = "a single numeric value",
    "integer" = "a single whole number",
    "bool" = "either TRUE or FALSE"
  )
  if (!check(x)) {
    cli::cli_abort(
      c(
        "i" = "{.arg {arg}} must be {msg}.",
        "x" = ifelse(
          length(x) == 1,
          "You've supplied a {.cls {class(x)}}.",
          "You've supplied a {.cls {class(x)}} of length {length(x)}."
        )
      ),
      call = call
    )
  }
}
#' @rdname type_check_helpers
#'
assert_vector <- function(
  x,
  type = c("character", "numeric", "integer", "logical"),
  n = 2,
  allow_null = FALSE,
  arg = caller_arg(x),
  call = caller_env()
) {
  if (allow_null && is.null(x)) {
    return(invisible(NULL))
  }
  type <- match.arg(type, choices = c("character", "numeric", "integer", "logical"))
  check <- switch(type,
    "character" = \(x) {
      is.character(x)
    },
    "numeric" = \(x) {
      is.numeric(x)
    },
    "integer" = \(x) {
      is_vector_wholenumber(x)
    },
    "logical" = \(x) {
      is.logical(x)
    }
  )
  msg <- switch(type,
    "character" = "a {.cls character} vector",
    "numeric" = "a {.cls numeric} vector",
    "integer" = "an {.cls integer} vector",
    "logical" = "a {.cls logical} vector"
  )
  if (!check(x) || length(x) < n) {
    cli::cli_abort(
      c(
        "i" = "{.arg {arg}} must be a {.cls {type}} vector with at least {n} element(s)",
        "x" = ifelse(
          length(x) >= n,
          "You've supplied a {.cls {class(x)}}.",
          "You've supplied a {.cls {class(x)}} of length {length(x)}."
        )
      ),
      call = call
    )
  }
}
#' @rdname type_check_helpers
#'
assert_class <- function(
  x,
  classes,
  allow_null = FALSE,
  arg = caller_arg(x),
  call = caller_env()
) {
  if (allow_null && is.null(x)) {
    return(invisible(NULL))
  }
  if (!inherits(x, classes)) {
    cli::cli_abort(
      c(
        "i" = "{.arg {arg}} must be a {.cls {classes}} object.",
        "x" = "You've supplied a {.cls {class(x)}} object."
      ),
      call = call
    )
  }
}
#' @rdname type_check_helpers
#'
assert_mpx_assay <- function(
  x,
  allow_null = FALSE,
  arg = caller_arg(x),
  call = caller_env()
) {
  if (allow_null && is.null(x)) {
    return(invisible(NULL))
  }
  if (!inherits(x, "MPXAssay")) {
    cli::cli_abort(
      c(
        "i" = "The selected Assay must be a {.cls {c('CellGraphAssay', 'CellGraphAssay5')}} object",
        "x" = "Got a {.cls {class(x)}} object."
      )
    )
  }
}
#' @rdname type_check_helpers
#'
assert_within_limits <- function(
  x,
  limits,
  arg = caller_arg(x),
  call = caller_env()
) {
  stopifnot(
    "`x` must be a single numeric value" =
      is.numeric(x),
    "`limits` must be a numeric vector of length 2" =
      is.numeric(limits) && length(limits) == 2
  )
  values_in_range <- dplyr::between(x, limits[1], limits[2])
  if (!all(values_in_range)) {
    cli::cli_abort(
      c(
        "i" = "Values in {.arg {arg}} must lie between {.val {limits[1]}} and {.val {limits[2]}}.",
        "x" = "{.val {x[!values_in_range]}} is outside of the range."
      ),
      call = call
    )
  }
}
#' @rdname type_check_helpers
#'
assert_function <- function(
  x,
  allow_null = FALSE,
  arg = caller_arg(x),
  call = caller_env()
) {
  if (allow_null && is.null(x)) {
    return(invisible(NULL))
  }
  if (!rlang::is_function(x)) {
    cli::cli_abort(
      c(
        "i" = "{.arg {arg}} must be a {.cls function}.",
        "x" = "You've supplied a {.cls {class(x)}}."
      ),
      call = call
    )
  }
}
#' @rdname type_check_helpers
#'
assert_file_exists <- function(
  x,
  allow_null = FALSE,
  call = caller_env()
) {
  if (allow_null && is.null(x)) {
    return(invisible(NULL))
  }
  if (!fs::file_exists(x)) {
    cli::cli_abort(
      c("x" = "File {.file {x}} doesn't exist."),
      call = call
    )
  }
}
#' @rdname type_check_helpers
#'
assert_file_ext <- function(
  x,
  ext,
  allow_null = FALSE,
  call = caller_env()
) {
  if (allow_null && is.null(x)) {
    return(invisible(NULL))
  }
  stopifnot(
    "`ext` must be a string" =
      rlang::is_string(ext)
  )
  if (fs::path_ext(x) != ext) {
    cli::cli_abort(
      c("x" = "File {.file {x}} is not a PXL file"),
      call = call
    )
  }
}
#' @rdname type_check_helpers
#'
assert_x_in_y <- function(
  x,
  y,
  allow_null = FALSE,
  arg_x = caller_arg(x),
  arg_y = caller_arg(y),
  call = caller_env()
) {
  if (allow_null && is.null(x)) {
    return(invisible(NULL))
  }
  stopifnot(
    "`x` and `y` must be vectors.`" =
      is_vector(x) && is_vector(y)
  )
  if (!all(x %in% y)) {
    missing_elements <- match(x, y)
    missing_elements <- which(is.na(missing_elements))

    cli::cli_abort(
      c(
        "i" = "All elements of {.arg {arg_x}} must be present in {.arg {arg_y}}.",
        "x" = "Positions of missing element(s) in {.arg {arg_x}}: {.val {missing_elements}}"
      ),
      call = call
    )
  }
}
#' @rdname type_check_helpers
#'
assert_single_values_are_different <- function(
  x,
  y,
  allow_null = FALSE,
  arg_x = caller_arg(x),
  arg_y = caller_arg(y),
  call = caller_env()
) {
  if (allow_null && is.null(x)) {
    return(invisible(NULL))
  }
  stopifnot(
    "`x` and `y` must be vectors" =
      is_vector(x) && is_vector(y),
    "`x` and `y` must have 1 element" =
      length(x) == 1 && length(y) == 1
  )
  if (x == y) {
    cli::cli_abort(
      c(
        "i" = "{.arg {arg_x}} and {.arg {arg_y}} must be different.",
        "x" = "{.arg {arg_x}={x}} and {.arg {arg_y}={y}}."
      ),
      call = call
    )
  }
}
#' @rdname type_check_helpers
#'
assert_singles_match <- function(
  x,
  y,
  arg_x = caller_arg(x),
  arg_y = caller_arg(y),
  call = caller_env()
) {
  stopifnot(
    "`x` and `y` must be numeric values`" =
      is.numeric(x) && (length(x) == 1) &&
        is.numeric(y) && (length(y) == 1)
  )
  if (x != y) {
    cli::cli_abort(
      c(
        "i" = "{.arg {arg_x}} and {.arg {arg_y}} must be identical.",
        "x" = "{.arg {arg_x}={x}}",
        "x" = "{.arg {arg_y}={y}}"
      ),
      call = call
    )
  }
}
#' @rdname type_check_helpers
#'
assert_length <- function(
  x,
  n = 1,
  allow_null = FALSE,
  arg_x = caller_arg(x),
  call = caller_env()
) {
  if (allow_null && is.null(x)) {
    return(invisible(NULL))
  }
  stopifnot(
    "`x` must be a vector" =
      is_vector(x),
    "`n` must be a numeric value" =
      (is_double(n) || is_integer(n))
  )
  if (length(x) != n) {
    cli::cli_abort(
      c("x" = "{.arg {arg_x}} must have exactly {.val {n}} element(s) but has {.val {length(x)}}."),
      call = call
    )
  }
}
#' @rdname type_check_helpers
#'
assert_max_length <- function(
  x,
  n = 1,
  allow_null = FALSE,
  arg_x = caller_arg(x),
  call = caller_env()
) {
  if (allow_null && is.null(x)) {
    return(invisible(NULL))
  }
  stopifnot(
    "`x` must be a vector" =
      is_vector(x),
    "`n` must be a numeric value" =
      (is_double(n) || is_integer(n))
  )
  if (length(x) > n) {
    cli::cli_abort(
      c(
        "i" = "{.arg x} cannot have more than {.val {n}} element(s)",
        "x" = "You've provided a {.cls {class(x)}} with {.val {length(x)}} elements"
      ),
      call = call
    )
  }
}
#' @rdname type_check_helpers
#'
assert_vectors_x_y_length_equal <- function(
  x,
  y,
  arg_x = caller_arg(x),
  arg_y = caller_arg(y),
  call = caller_env()
) {
  stopifnot(
    "`x` and `y` must be vectors.`" =
      is_vector(x) && is_vector(y)
  )
  if (length(x) != length(y)) {
    cli::cli_abort(
      c(
        "i" = "{.arg {arg_x}} and {.arg {arg_y}} must have the same lengths.",
        "x" = "Length of {.arg {arg_x}}: {length(x)}",
        "x" = "Length of {.arg {arg_y}}: {length(y)}"
      ),
      call = call
    )
  }
}
#' @rdname type_check_helpers
#'
assert_unique <- function(
  x,
  arg = caller_arg(x),
  call = caller_env()
) {
  stopifnot(
    "`x` must be a vector.`" =
      is_vector(x)
  )
  if (length(x) != length(unique(x))) {
    dup_elements <- which(duplicated(x))
    cli::cli_abort(
      c(
        "i" = "{.arg {arg}} must have unique elements.",
        "x" = "{.arg {arg}} has duplicated elements at positions: {.val {dup_elements}}."
      ),
      call = call
    )
  }
}
#' @rdname type_check_helpers
#'
assert_col_class <- function(
  x,
  data,
  classes,
  allow_null = FALSE,
  arg_data = caller_arg(data),
  call = caller_env()
) {
  if (allow_null && is.null(x)) {
    return(invisible(NULL))
  }
  stopifnot(
    "`x` must be a string`" =
      is_string(x),
    "`data` must be a data.frame-like object`" =
      inherits(data, "data.frame")
  )
  col_x <- data[, x, drop = TRUE]
  if (!inherits(col_x, classes)) {
    cli::cli_abort(
      c(
        "i" = ifelse(
          length(classes) == 1,
          "Column {.val {x}} in {.arg {arg_data}} must be a {.cls {classes}}.",
          "Column {.val {x}} in {.arg {arg_data}} must be one of {.cls {classes}}."
        ),
        "x" = "Column {.val {x}} is a {.cls {class(col_x)}}."
      ),
      call = call
    )
  }
}
#' @rdname type_check_helpers
#'
assert_col_in_data <- function(
  x,
  data,
  allow_null = FALSE,
  arg_data = caller_arg(data),
  call = caller_env()
) {
  if (allow_null && is.null(x)) {
    return(invisible(NULL))
  }
  stopifnot(
    "`x` must be a string`" =
      is_string(x),
    "`data` must be a data.frame-like object`" =
      inherits(data, "data.frame")
  )
  if (!(x %in% names(data))) {
    cli::cli_abort(
      c("x" = "Column {.str {x}} is missing from {.arg {arg_data}}."),
      call = call
    )
  }
}
#' @rdname type_check_helpers
#'
assert_non_empty_object <- function(
  x,
  classes,
  allow_null = FALSE,
  arg = caller_arg(x),
  call = caller_env()
) {
  if (allow_null) {
    return(invisible(NULL))
  }
  if (!(inherits(x, what = classes) && length(x) > 0)) {
    cli::cli_abort(
      c(
        "i" = "{.arg {arg}} must be a non-empty {.cls {classes}} object.",
        "x" = ifelse(
          length(x) == 0,
          "You've supplied an empty {.cls {class(x)}}.",
          "You've supplied a {.cls {class(x)}}."
        )
      ),
      call = call
    )
  }
}
#' @rdname type_check_helpers
#'
assert_vectors_match <- function(
  x,
  y,
  allow_null = FALSE,
  arg_x = caller_arg(x),
  arg_y = caller_arg(y),
  call = caller_env()
) {
  if (allow_null) {
    return(invisible(NULL))
  }
  stopifnot(
    "`x` and `y` must be vectors.`" =
      is_vector(x) && is_vector(y),
    "`x` and `y` cannot have NA values" =
      !anyNA(x) && !anyNA(y)
  )
  if (length(x) != length(y)) {
    cli::cli_abort(
      c(
        "i" = "{.arg {arg_x}} and {.arg {arg_y}} must have the same lengths.",
        "x" = "Length of {.arg {arg_x}}: {length(x)}",
        "x" = "Length of {.arg {arg_y}}: {length(y)}"
      ),
      call = call
    )
  }
  if (!all(x == y)) {
    unmatched_elements <- which(x != y)
    cli::cli_abort(
      c(
        "i" = "{.arg {arg_x}} and {.arg {arg_y}} must be identical.",
        "x" = "{.arg {arg_x}} and {.arg {arg_y}} differ at positions: {.val {unmatched_elements}}."
      ),
      call = call
    )
  }
}


is_scalar_wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  if (!length(x) == 1) {
    return(FALSE)
  }
  if (rlang::is_double(x) || rlang::is_integer(x)) {
    abs(x - round(x)) < tol
  } else {
    FALSE
  }
}

is_vector_wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  if (rlang::is_double(x) || rlang::is_integer(x)) {
    all(abs(x - round(x)) < tol)
  } else {
    FALSE
  }
}
