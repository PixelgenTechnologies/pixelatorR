pxl_file <- minimal_pna_pxl_file()

test_that("ReadPNA_edgelist works as expected", {
  # tbl_df
  expect_no_error(el <- ReadPNA_edgelist(pxl_file, lazy = FALSE, umi_data_type = "string"))
  expected_data <-
    structure(
      list(
        marker_1 = c("CD6", "CD6"),
        marker_2 = c("CD44",
                     "CD44"),
        umi1 = c("10004431758516698", "10004431758516698"),
        umi2 = c("10532227491147037", "66853920218164601"),
        read_count = 2:1,
        uei_count = c(1L, 1L),
        component = c("c3c393e9a17c1981",
                      "c3c393e9a17c1981")
      ),
      row.names = c(NA,-2L),
      class = c("tbl_df",
                "tbl", "data.frame")
    )
  expect_identical(head(el %>% arrange(umi1), 2), expected_data)

  # tbl_lazy
  expect_no_error(el <- ReadPNA_edgelist(pxl_file))
  expect_equal(c("marker_1", "marker_2", "umi1", "umi2", "read_count", "uei_count", "component"), colnames(el))
  expect_s3_class(el, "tbl_lazy")
  expect_equal(dim(el %>% collect()), c(528594, 7))

  expect_no_error(el <- ReadPNA_edgelist(pxl_file, cells = "c3c393e9a17c1981"))
  expect_equal(dim(el %>% collect()), c(110657, 7))
})

test_that("ReadPNA_edgelist fails with invalid input", {
  expect_error(ReadPNA_edgelist("Invalid"))
  expect_error(ReadPNA_edgelist(pxl_file, cells = "Invalid"))
  expect_error(ReadPNA_edgelist(pxl_file, umi_data_type = "Invalid"))
  expect_error(ReadPNA_edgelist(pxl_file, verbose = "Invalid"))
  expect_error(ReadPNA_edgelist(pxl_file, lazy = "Invalid"))
})
