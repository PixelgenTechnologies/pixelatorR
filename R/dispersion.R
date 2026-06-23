#' @rdname CalculateDispersion
#' @method CalculateDispersion matrix
#'
#' @export
#'
CalculateDispersion.matrix <- function(
  object,
  method = c("gini", "tau"),
  margin = 2,
  ...
) {
  # Validate input parameters
  assert_class(object, "matrix")
  assert_non_empty_object(object, "matrix")
  assert_vector(method, type = "character", n = 1)
  assert_single_value(margin, type = "numeric")
  assert_x_in_y(method, c("gini", "tau"))
  assert_x_in_y(margin, c(1, 2))

  method <- match.arg(method, c("gini", "tau"))

  dispersion_function <-
    switch(method,
      "gini" = function(x) {
        n <- length(x)
        if (n <= 1L) {
          return(0)
        }
        sum_x <- sum(x)
        if (sum_x == 0) {
          return(0)
        }
        x_sorted <- sort(x)
        i <- seq_len(n)
        (2 * sum(i * x_sorted) / (n * sum_x)) - (n + 1) / n
      },
      "tau" = function(x) {
        n <- length(x)
        if (n <= 1L) {
          return(0)
        }
        max_x <- max(x)
        if (max_x == 0) {
          return(0)
        }
        sum(1 - x / max_x) / (n - 1)
      }
    )

  object %>%
    apply(margin, dispersion_function)
}

#' @rdname CalculateDispersion
#' @method CalculateDispersion data.frame
#'
#' @export
#'
CalculateDispersion.data.frame <- function(
  object,
  method = c("gini", "tau"),
  margin = 2,
  ...
) {
  for (colmn in colnames(object)) {
    assert_col_class(colmn, object, "numeric")
  }

  CalculateDispersion(
    object = as.matrix(object),
    method = method,
    margin = margin,
    ...
  )
}

#' @rdname CalculateDispersion
#' @method CalculateDispersion Matrix
#'
#' @export
#'
CalculateDispersion.Matrix <- function(
  object,
  method = c("gini", "tau"),
  margin = 2,
  ...
) {
  CalculateDispersion(
    object = as.matrix(object),
    method = method,
    margin = margin,
    ...
  )
}

#' @rdname CalculateDispersion
#' @method CalculateDispersion Seurat
#'
#' @export
#'
CalculateDispersion.Seurat <- function(
  object,
  method = c("gini", "tau"),
  margin = 2,
  assay = NULL,
  layer = NULL,
  metadata_name = NULL,
  ...
) {
  assert_single_value(metadata_name, type = "string", allow_null = TRUE)

  method <- match.arg(method, c("gini", "tau"))

  dispersion <-
    CalculateDispersion(
      object = LayerData(object,
        assay = assay,
        layer = layer
      ),
      method = method,
      margin = margin,
      ...
    )

  if (margin == 2) {
    if (is.null(metadata_name)) {
      metadata_name <- paste0("dispersion_", method)
    }
    return(AddMetaData(object, dispersion, col.name = metadata_name))
  } else {
    return(dispersion)
  }
}
