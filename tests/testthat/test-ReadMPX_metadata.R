pxl_file <- system.file("extdata/five_cells",
                        "five_cells.pxl",
                        package = "pixelatorR")
meta_data <- ReadMPX_metadata(pxl_file)

test_that("ReadMPX_metadata works as expected", {

  expected_data <-
    structure(
      list(
        version = "0.12.0",
        sample = "Mock_data",
        analysis = list(
          params = c(
            compute_polarization = "TRUE",
            compute_colocalization = "TRUE",
            use_full_bipartite = "FALSE",
            polarization_normalization = "clr",
            polarization_binarization = "FALSE",
            colocalization_transformation = "log1p",
            colocalization_neighbourhood_size = "1",
            colocalization_n_permutations = "50",
            colocalization_min_region_count = "5"
          )
        )
      ),
      row.names = c(NA,-1L),
      class = c("pixelator_metadata", "tbl_df", "tbl", "data.frame")
    )
  expect_equal(meta_data, expected_data)

  # print method
  captured_output <- capture_output_lines(print(meta_data))
  expect_equal(
    captured_output,
    c(
      "Sample 1 name: \t\tMock_data",
      "Pixelator version: \t0.12.0",
      "",
      "# A tibble: 9 × 2",
      "  parameter                         value",
      "  <chr>                             <chr>",
      "1 compute_polarization              TRUE ",
      "2 compute_colocalization            TRUE ",
      "3 use_full_bipartite                FALSE",
      "4 polarization_normalization        clr  ",
      "5 polarization_binarization         FALSE",
      "6 colocalization_transformation     log1p",
      "7 colocalization_neighbourhood_size 1    ",
      "8 colocalization_n_permutations     50   ",
      "9 colocalization_min_region_count   5    "
    )
  )

})
