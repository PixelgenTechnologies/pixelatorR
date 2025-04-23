pxl_file <- minimal_pna_pxl_file()

test_that("ReadPNA_counts works as expected", {
  # matrix
  expect_no_error(counts <- ReadPNA_counts(pxl_file))
  expected_data <-
    new(
      "dgCMatrix",
      i = c(0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L),
      p = c(0L, 2L, 4L, 6L, 8L, 10L),
      Dim = c(2L, 5L),
      Dimnames = list(
        c("HLA-ABC", "B2M"),
        c(
          "0a45497c6bfbfb22",
          "2708240b908e2eba",
          "c3c393e9a17c1981",
          "d4074c845bb62800",
          "efe0ed189cb499fc"
        )
      ),
      x = c(
        865, 1182, 2077, 3448, 2480, 5307, 2994, 9753,
        2212, 6082
      ),
      factors = list()
    )
  expect_identical(head(counts, 2), expected_data)
})

test_that("ReadPNA_counts fails with invalid input", {
  expect_error(ReadPNA_counts("Invalid"))
  expect_error(ReadPNA_counts(pxl_file, verbose = "Invalid"))
  expect_error(ReadPNA_counts(pxl_file, return_list = "Invalid"))
})
