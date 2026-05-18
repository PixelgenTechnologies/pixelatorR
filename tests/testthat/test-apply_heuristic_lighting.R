test_that("apply_heuristic_lighting works as expected", {
  set.seed(123)

  layout <- tibble::tibble(
    x = rnorm(100),
    y = rnorm(100),
    z = rnorm(100)
  )

  expect_no_error(
    illum <- apply_heuristic_lighting(layout)
  )

  expect_equal(sum(illum), 56.83794, tolerance = 1e-4)
})
