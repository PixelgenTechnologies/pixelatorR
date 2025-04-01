pxl_file <- system.file("extdata/five_cells",
  "five_cells.pxl",
  package = "pixelatorR"
)
se <- ReadMPX_Seurat(pxl_file)

test_that("NormalizeMethods work as expected", {
  expect_no_error(
    .normalize_method_dsb(LayerData(se, "counts"),
      isotype_controls = c("mIgG1", "mIgG2a", "mIgG2b")
    )
  )
  expect_no_error(
    .normalize_method_clr(LayerData(se, "counts"))
  )
})

test_that("Normalize works as expected", {
  # Valid input
  expect_no_error(
    Normalize(LayerData(se, "counts"))
  )
  expect_no_error(
    Normalize(se[["mpxCells"]])
  )
  expect_no_error(
    norm_data_clr <- Normalize(se, method = "clr")
  )

  expect_equal(
    head(LayerData(norm_data_clr)),
    new("dgCMatrix",
      i = c(
        0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L,
        4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L,
        2L, 3L, 4L, 5L
      ), p = c(0L, 6L, 12L, 18L, 24L, 30L), Dim = 6:5,
      Dimnames = list(c(
        "CD274", "CD44", "CD25", "CD279", "CD41",
        "HLA-ABC"
      ), c(
        "RCVCMP0000217", "RCVCMP0000118", "RCVCMP0000487",
        "RCVCMP0000655", "RCVCMP0000263"
      )), x = c(
        -0.0899921902406002,
        1.73625345505862, -0.32638096830483, 0.654448284706896, -1.93581888073893,
        5.52082822281565, -0.728278799753101, 1.23783405661973, -0.88242947958036,
        -1.98104176824847, 0.270250030358026, 5.90522758578796, -0.799819236234203,
        1.51271618761301, -0.46334699961299, -1.02296278754841, -0.157965350061808,
        6.04903493482777, 0.509830230167604, 3.15777650720011, -0.365638507186296,
        -0.470999022844122, -0.365638507186296, 5.51461511058568,
        -0.128021925145374, 2.77196529145023, 0.291831920414889,
        -2.16490385240641, -2.16490385240641, 6.04422313688879
      ),
      factors = list()
    )
  )

  expect_no_error(
    norm_data_dsb <- Normalize(se)
  )

  expect_equal(
    head(LayerData(norm_data_dsb)),
    new("dgCMatrix",
      i = c(
        0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L,
        4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L,
        2L, 3L, 4L, 5L
      ), p = c(0L, 6L, 12L, 18L, 24L, 30L), Dim = 6:5,
      Dimnames = list(c(
        "CD274", "CD44", "CD25", "CD279", "CD41",
        "HLA-ABC"
      ), c(
        "RCVCMP0000217", "RCVCMP0000118", "RCVCMP0000487",
        "RCVCMP0000655", "RCVCMP0000263"
      )), x = c(
        0.583158270809046,
        0.079078082542696, -0.154725746414455, 1.49282321117081,
        0.614776868067649, -0.136067006221542, 0.141674339857816,
        -0.171753260080644, 0.111411814794127, 0.101932380510186,
        1.01582427804594, 0.480672998470627, 0.059390105332717, 0.0523654846610522,
        0.0269396277277562, 0.223589218874532, 2.15433057618296,
        0.585732137010824, 1.23777961281461, 1.59378291104827, 0.333471374630927,
        1.21408432235121, 0.726806355392755, -0.0606221174819132,
        0.863241929077271, 1.44101996837564, 0.881550999939355, -0.840953856900626,
        0.383868338140238, 0.711121074357595
      ), factors = list()
    )
  )

  # Invalid input
  expect_error(
    Normalize(se,
      method = "invalid_method"
    ),
    regexp = "'arg' should be"
  )
  expect_error(
    Normalize(se,
      method = "dsb",
      isotype_controls = c(
        "isotype_control_that_definitely_does_not_exist",
        "isotype_control_that_could_exist_but_does_not"
      )
    ),
    regexp = "All elements of `isotype_controls` must be present in"
  )

  expect_error(
    Normalize(se,
      method = "dsb",
      isotype_controls = NULL
    )
  )
})

# Test Assay5 only
se5 <- se
se <- Normalize(se, method = "clr")
suppressWarnings({
  se5[["mpxCells"]] <- as(object = se[["mpxCells"]], Class = "Assay5")
})
se5 <- Normalize(se5, method = "clr")
test_that("CellGraphAssay5 and Assay5 are equal", {
  expect_equal(
    se[["mpxCells"]]$data,
    se5[["mpxCells"]]$data
  )
})
