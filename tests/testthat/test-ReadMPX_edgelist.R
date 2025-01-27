pxl_file <- system.file("extdata/five_cells", "five_cells.pxl", package = "pixelatorR")

test_that("ReadMPX_arrow_edgelist works as expected", {
  el <- ReadMPX_arrow_edgelist(pxl_file)
  expect_type(el, type = "environment")
  expect_equal(names(el), c(
    "upia", "upib", "marker", "count",
    "component"
  ))
  expect_no_error({
    msg <- capture_messages({
      el <- ReadMPX_arrow_edgelist(pxl_file)
    })
  })

  tmp_f <- fs::path_temp("test.parquet")
  expect_no_error({
    el <-
      ReadMPX_arrow_edgelist(
        pxl_file,
        edge_list_file = tmp_f,
        verbose = FALSE
      )
  })
  expect_true(fs::file_exists(tmp_f))
  expect_true(basename(tmp_f) == "test.parquet")
})


test_that("ReadMPX_arrow_edgelist fails when invalid input is provided", {
  expect_error({
    el <- ReadMPX_arrow_edgelist("Invalid file",
      outdir = tempdir()
    )
  })
  expect_error({
    el <- ReadMPX_arrow_edgelist("Invalid file", edge_list_file = fs::path_temp("test.pxl")
    )
  })
})
