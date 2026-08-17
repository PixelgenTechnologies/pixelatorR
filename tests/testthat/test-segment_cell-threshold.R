test_that(".kmeans_midpoint returns center midpoint", {
  expect_equal(pixelatorR:::.kmeans_midpoint(c(0.2, 0.8)), 0.5)
  expect_equal(pixelatorR:::.kmeans_midpoint(c(0.65, 0.35)), 0.5)
  expect_equal(pixelatorR:::.kmeans_midpoint(c(0.1, 0.3)), 0.2)
})

test_that(".kmeans_midpoint rejects invalid center vectors", {
  expect_error(pixelatorR:::.kmeans_midpoint(0.5))
  expect_error(pixelatorR:::.kmeans_midpoint(c("a", "b")))
})
