pxl_file <- minimal_pna_pxl_file()

test_that("ReadPNA_edgelist works as expected", {
  # tbl_df
  expect_no_error(el <- ReadPNA_edgelist(pxl_file, umi_data_type = "string"))
  expected_data <-
    structure(
      list(
        umi1 = c("10004431758516698", "10004431758516698"),
        umi2 = c("10532227491147037", "66853920218164601"),
        marker_1 = c("CD6",
                     "CD6"),
        marker_2 = c("CD44", "CD44"),
        component = c("c3c393e9a17c1981",
                      "c3c393e9a17c1981"),
        read_count = 2:1,
        uei_count = c(1L, 1L)
      ),
      row.names = c(NA,-2L),
      class = c("tbl_df", "tbl", "data.frame")
    )
  expect_identical(head(el %>% arrange(umi1), 2), expected_data)

  # Table
  expect_no_error(el <- ReadPNA_edgelist(pxl_file))
  expect_equal(dim(el), c(528594, 7))
  expect_s3_class(el, "tbl_df")
})

test_that("ReadPNA_edgelist fails with invalid input", {
  expect_error(ReadPNA_edgelist("Invalid"))
  expect_error(ReadPNA_edgelist(pxl_file, verbose = "Invalid"))
  expect_error(ReadPNA_edgelist(pxl_file, return_tibble = "Invalid"))
})
