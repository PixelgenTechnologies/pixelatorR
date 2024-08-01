#' The pixelatorR style
#'
#' Style code according to the pixelatorR style guide. This style guide is
#' based on the tidyverse style guide, but with some modifications.
#'
#' - \code{unindent_fun_dec} is set to \code{NULL}
#' - \code{remove_line_breaks_in_fun_dec} is set to \code{NULL}
#'
#' @return A list of transformers that can be used with e.g. \code{styler::style_text()}
#'
#' @export
#'
pixelatorR_style <- function() {
  expect_styler()
  transformers <- styler::tidyverse_style()
  transformers$indention$unindent_fun_dec <- NULL
  transformers$line_break$remove_line_breaks_in_fun_dec <- NULL
  return(transformers)
}
