library(dplyr)
library(tibble)

# Generate random points that are offset to (100, 100, 100)
xyz <- matrix(rnorm(600, mean = 100, sd = 20), ncol = 3,
              dimnames = list(NULL, c("x", "y", "z"))) %>%
  as_tibble()

test_that("Layout coordinate utility functions work as expected", {

  # project_layout_coordinates_on_unit_sphere
  xyz %>% project_layout_coordinates_on_unit_sphere() %>%
    mutate(across(x:z, ~.x^2)) %>%
    rowSums() %>%
    sqrt() %>%
    expect_equal(., rep(1, nrow(xyz)))

  # normalize_layout_coordinates
  xyz %>% normalize_layout_coordinates() %>%
    mutate(across(x:z, ~.x^2)) %>%
    rowSums() %>%
    sqrt() %>%
    median() %>%
    expect_equal(., 1)

  # center_layout_coordinates
  xyz %>% center_layout_coordinates() %>%
    apply(2, mean) %>%
    near(., y = 0, tol = 1e-10) %>%
    all() %>%
    expect_true()
})

test_that("Layout coordinate utility functions fails with invalid input", {

  # project_layout_coordinates_on_unit_sphere
  expect_error(project_layout_coordinates_on_unit_sphere("Invalid"),
               "'layout' must be a non-empty, matrix-like object")
  expect_error(project_layout_coordinates_on_unit_sphere(xyz[, 1:2]),
               "'layout' can only have 3 columns")

  # normalize_layout_coordinates
  expect_error(normalize_layout_coordinates("Invalid"),
               "'layout' must be a non-empty, matrix-like object")

  # center_layout_coordinates
  expect_error(center_layout_coordinates("Invalid"),
               "'layout' must be a non-empty, matrix-like object")
})
