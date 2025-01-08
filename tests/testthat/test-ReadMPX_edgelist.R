test_that("ReadMPX_arrow_edgelist works as expected", {
  el <- ReadMPX_arrow_edgelist(system.file("extdata/five_cells", "five_cells.pxl", package = "pixelatorR"))
  expect_type(el, type = "environment")
  expect_equal(names(el), c(
    "upia", "upib", "marker", "count",
    "component"
  ))
  expect_no_error({
    msg <- capture_messages({
      el <- ReadMPX_arrow_edgelist(system.file("extdata/five_cells", "five_cells.pxl", package = "pixelatorR"))
    })
  })
})


test_that("ReadMPX_arrow_edgelist fails when invalid input is provided", {
  expect_error(
    {
      el <- ReadMPX_arrow_edgelist("Invalid file",
        outdir = tempdir()
      )
    }
  )
})
