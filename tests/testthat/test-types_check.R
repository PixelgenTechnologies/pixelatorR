test_that("assert_col_class works as expected", {
  # A column named "x" must not shadow the column selected with the `x` argument
  data <- tibble::tibble(
    x = c("a", "b"),
    log2_ratio = c(1.5, 2.5)
  )

  expect_no_error(assert_col_class("log2_ratio", data, classes = "numeric"))
  expect_no_error(assert_col_class("x", data, classes = c("character", "factor")))
  expect_no_error(assert_col_class(NULL, data, classes = "numeric", allow_null = TRUE))

  expect_error(assert_col_class("log2_ratio", data, classes = "character"), "must be a")
  expect_error(assert_col_class("x", data, classes = "numeric"), "must be a")
  expect_error(assert_col_class("missing_column", data, classes = "numeric"))
})
