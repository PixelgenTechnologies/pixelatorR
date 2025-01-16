options(Seurat.object.assay.version = "v3")

pxl_file <- system.file("extdata/five_cells",
  "five_cells.pxl",
  package = "pixelatorR"
)
# Seurat objects
se <- ReadMPX_Seurat(pxl_file)
se <- merge(se, rep(list(se), 9), add.cell.ids = LETTERS[1:10])
se$sample <- c("T", "R", "C", "C", "C") %>% rep(times = 10)
se <- Seurat::NormalizeData(se, normalization.method = "CLR", margin = 2)

test_that("RunDAA works as expected on a Seurat object", {
  # One target
  expect_no_error(suppressWarnings(daa_markers <- RunDAA(se,
    contrast_column = "sample",
    targets = "T", reference = "C"
  )))

  expected_result <-
    structure(
      list(
        marker = c("CD137", "CD62P"),
        p = c(
          3.73436610215361e-07,
          1.42675642066992e-06
        ),
        p_adj = c(2.98749288172289e-05, 0.000114140513653594),
        difference = c(0.0672793901912557, -0.130339927355884),
        pct_1 = c(
          1,
          0
        ),
        pct_2 = c(0.333, 1),
        target = c("T", "T"),
        reference = c(
          "C",
          "C"
        )
      ),
      row.names = c(NA, -2L),
      class = c("tbl_df", "tbl", "data.frame")
    )

  expect_equal(daa_markers[1:2, ], expected_result)

  # Multiple targets
  expect_no_error(suppressWarnings(daa_markers <- RunDAA(se,
    contrast_column = "sample",
    targets = c("T", "R"), reference = "C"
  )))

  expected_result <-
    structure(
      list(
        marker = c("CD137", "CD62P", "CD55", "CD314"),
        p = c(
          3.73436610215361e-07,
          1.42675642066992e-06,
          0.11041579736434,
          0.11041579736434
        ),
        p_adj = c(
          5.97498576344578e-05, 0.000228281027307187,
          1, 1
        ),
        difference = c(
          0.0672793901912557,
          -0.130339927355884,
          0.372004130942341,
          -0.345923096372618
        ),
        pct_1 = c(
          1, 0, 1,
          1
        ),
        pct_2 = c(0.333, 1, 1, 1),
        target = c(
          "T", "T", "R",
          "R"
        ),
        reference = c("C", "C", "C", "C")
      ),
      row.names = c(NA, -4L),
      class = c("tbl_df", "tbl", "data.frame")
    )

  expect_equal(daa_markers[c(1:2, (nrow(daa_markers) - 1):nrow(daa_markers)), ], expected_result)
})

test_that("RunDAA fails with invalid input", {
  expect_error(
    dpa_markers <- RunDAA(se),
    'argument "contrast_column" is missing, with no default'
  )
  expect_error(
    dpa_markers <- RunDAA(se, reference = "Sample1"),
    'argument "contrast_column" is missing, with no default'
  )
  expect_error(
    dpa_markers <- RunDAA(se, contrast_column = "sample", targets = "Sample1"),
    'argument "reference" is missing, with no default'
  )
  expect_error(
    dpa_markers <- RunDAA(se, contrast_column = "Invalid")
  )
  expect_error(
    dpa_markers <- RunDAA(se, contrast_column = "sample", targets = "Invalid")
  )
  expect_error(
    dpa_markers <- RunDAA(se, contrast_column = "sample", targets = "T", reference = "Invalid")
  )
  expect_error(
    dpa_markers <- RunDAA(se, contrast_column = "sample", targets = c("T", "C"), reference = "C")
  )
  expect_error(
    dpa_markers <- RunDAA(se, contrast_column = "sample", targets = "T", reference = "C", group_vars = "sample")
  )

  expect_error(
    dpa_markers <- RunDAA(se, contrast_column = "sample", targets = "T", reference = "C", group_vars = "Invalid")
  )
})
