# Type check helpers

Utility functions for type checking in `pixelatorR`. These functions
throw informative error messages if the checks fail.

## Usage

``` r
assert_single_value(
  x,
  type = c("string", "numeric", "integer", "bool"),
  allow_null = FALSE,
  arg = caller_arg(x),
  call = caller_env()
)

assert_vector(
  x,
  type = c("character", "numeric", "integer", "logical"),
  n = 2,
  allow_null = FALSE,
  arg = caller_arg(x),
  call = caller_env()
)

assert_class(
  x,
  classes,
  allow_null = FALSE,
  arg = caller_arg(x),
  call = caller_env()
)

assert_mpx_assay(
  x,
  allow_null = FALSE,
  arg = caller_arg(x),
  call = caller_env()
)

assert_pna_assay(
  x,
  allow_null = FALSE,
  arg = caller_arg(x),
  call = caller_env()
)

assert_pixel_assay(
  x,
  allow_null = FALSE,
  arg = caller_arg(x),
  call = caller_env()
)

assert_within_limits(x, limits, arg = caller_arg(x), call = caller_env())

assert_function(
  x,
  allow_null = FALSE,
  arg = caller_arg(x),
  call = caller_env()
)

assert_file_exists(x, allow_null = FALSE, call = caller_env())

assert_file_ext(x, ext, allow_null = FALSE, call = caller_env())

assert_x_in_y(
  x,
  y,
  allow_null = FALSE,
  arg_x = caller_arg(x),
  arg_y = caller_arg(y),
  call = caller_env()
)

assert_single_values_are_different(
  x,
  y,
  allow_null = FALSE,
  arg_x = caller_arg(x),
  arg_y = caller_arg(y),
  call = caller_env()
)

assert_singles_match(
  x,
  y,
  arg_x = caller_arg(x),
  arg_y = caller_arg(y),
  call = caller_env()
)

assert_length(
  x,
  n = 1,
  allow_null = FALSE,
  arg_x = caller_arg(x),
  call = caller_env()
)

assert_max_length(
  x,
  n = 1,
  allow_null = FALSE,
  arg_x = caller_arg(x),
  call = caller_env()
)

assert_vectors_x_y_length_equal(
  x,
  y,
  arg_x = caller_arg(x),
  arg_y = caller_arg(y),
  call = caller_env()
)

assert_unique(x, arg = caller_arg(x), call = caller_env())

assert_col_class(
  x,
  data,
  classes,
  allow_null = FALSE,
  arg_data = caller_arg(data),
  call = caller_env()
)

assert_col_in_data(
  x,
  data,
  allow_null = FALSE,
  arg_data = caller_arg(data),
  call = caller_env()
)

assert_non_empty_object(
  x,
  classes,
  allow_null = FALSE,
  arg = caller_arg(x),
  call = caller_env()
)

assert_is_one_of(
  x,
  choices,
  allow_null = FALSE,
  arg = caller_arg(x),
  call = caller_env()
)

assert_different(
  x,
  y,
  allow_null = FALSE,
  arg_x = caller_arg(x),
  arg_y = caller_arg(y),
  call = caller_env()
)

assert_vectors_match(
  x,
  y,
  allow_null = FALSE,
  arg_x = caller_arg(x),
  arg_y = caller_arg(y),
  call = caller_env()
)

assert_valid_color(
  x,
  n = 1,
  allow_null = FALSE,
  arg = caller_arg(x),
  call = caller_env()
)
```

## Arguments

- x, y:

  An object to check

- type:

  A character string specifying the type or vector class to check
  against.

- allow_null:

  Either `TRUE` or `FALSE`. If `TRUE`, `x` can be `NULL` and the check
  passes.

- arg, arg_x, arg_y, arg_data:

  The name of an argument to check. Used for error messages.

- call:

  An environment, typically the environment in which the function was
  called

- n:

  An integer

- classes:

  A character vector of classes to check against

- limits:

  A numeric vector of length 2 specifying a lower and upper limit

- ext:

  A character string specifying a file extension

- data:

  A data-frame like object

- choices:

  A vector of allowed values

## Value

Nothing if the check passes, otherwise throws an error.

## Details

- `assert_single_value` checks if `x` is a single value of a specified
  `type`.

- `assert_vector` checks if `x` is a vector of a specified `type` and
  with at least `n` elements.

- `assert_class` checks if `x` is of a specific class.

- `assert_mpx_assay` checks if `x` is a `CellGraphAssay` or
  `CellGraphAssay5`.

- `assert_pna_assay` checks if `x` is a `PNAAssay` or `PNAAssay5`.

- `assert_pixel_assay` checks if `x` is a `CellGraphAssay`,
  `CellGraphAssay5`, `PNAAssay` or `PNAAssay5`.

- `assert_x_in_y` checks if all elements of `x` are in `y`.

- `assert_single_values_are_different` checks if `x` and `y` are not
  different strings.

- `assert_col_class` checks if column `x` in `data` is of a specific
  class.

- `assert_col_in_data` checks if column `x` is present in `data`.

- `assert_vectors_match` checks if vectors `x` and `y` are identical.

- `assert_length` checks if vector `x` has a specific length.

- `assert_within_limits` checks if values in `x` are within `limits`.

- `assert_function` checks `x` is a function.

- `assert_file_exists` checks is file `x` exists.

- `assert_file_ext` checks is file `x` has file extension `ext`.

- `assert_unique` checks if values in `x` are unique.
