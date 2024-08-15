options(Seurat.object.assay.version = "v3")

pxl_file <- system.file("extdata/five_cells",
                        "five_cells.pxl",
                        package = "pixelatorR")

# Load polarization scores
polarization_table1 <- polarization_table2 <- ReadMPX_polarization(pxl_file)
polarization_table1$sample <- "Sample1"
polarization_table2$sample <- "Sample2"
polarization_table_merged <-  bind_rows(polarization_table1, polarization_table2)

# Seurat objects
seur1 <- seur2 <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)
seur1$sample <- "Sample1"
seur2$sample <- "Sample2"
seur_merged <- merge(seur1, seur2, add.cell.ids = c("Sample1", "Sample2"))

test_that("RunDPA works as expected on a data.frame", {

  # morans_z
  expect_no_error(suppressWarnings(dpa_markers <- RunDPA(polarization_table_merged,
    contrast_column = "sample",
    targets = "Sample1", reference = "Sample2"
  )))

  expected_result <-
    structure(
      list(
        estimate = c(
          `difference in location` = 0,
          `difference in location` = 0
        ),
        data_type = c("morans_z", "morans_z"),
        target = c("Sample1",
                   "Sample1"),
        reference = c("Sample2", "Sample2"),
        n1 = 4:5,
        n2 = 4:5,
        statistic = c(W = 8, W = 12.5),
        p = c(1, 1),
        p_adj = c(1,
                  1),
        conf.low = c(-0.159036147906079,-8.12133146083903),
        conf.high = c(0.159036147906079, 8.12133146083903),
        method = c("Wilcoxon",
                   "Wilcoxon"),
        alternative = c("two.sided", "two.sided"),
        marker = c("ACTB",
                   "B2M")
      ),
      row.names = c(NA,-2L),
      class = c("tbl_df", "tbl",
                "data.frame")
    )

  expect_equal(dpa_markers[1:2, ], expected_result)

  # morans_i
  expect_no_error(suppressWarnings(
    dpa_markers <- RunDPA(
      polarization_table_merged,
      contrast_column = "sample",
      targets = "Sample1",
      reference = "Sample2",
      polarity_metric = "morans_i"
    )
  ))

  expected_result <-
    structure(
      list(
        estimate = c(
          `difference in location` = 0,
          `difference in location` = 0
        ),
        data_type = c("morans_i", "morans_i"),
        target = c("Sample1",
                   "Sample1"),
        reference = c("Sample2", "Sample2"),
        n1 = 4:5,
        n2 = 4:5,
        statistic = c(W = 8, W = 12.5),
        p = c(1, 1),
        p_adj = c(1,
                  1),
        conf.low = c(-0.00179887506379854,-0.052500600183388),
        conf.high = c(0.00179887506379854, 0.052500600183388),
        method = c("Wilcoxon", "Wilcoxon"),
        alternative = c("two.sided",
                        "two.sided"),
        marker = c("ACTB", "B2M")
      ),
      row.names = c(NA,-2L),
      class = c("tbl_df", "tbl", "data.frame")
    )

  expect_equal(dpa_markers[1:2, ], expected_result)

  # Automatic selection of targets
  expect_no_error(suppressWarnings(dpa_markers <- RunDPA(seur_merged, contrast_column = "sample", reference = "Sample2")))

})

test_that("RunDPA works as expected on a Seurat object", {

  expect_no_error(suppressWarnings(dpa_markers <- RunDPA(seur_merged,
    contrast_column = "sample",
    targets = "Sample1", reference = "Sample2"
  )))

  expected_result <-
    structure(
      list(
        estimate = c(
          `difference in location` = 0,
          `difference in location` = 0
        ),
        data_type = c("morans_z", "morans_z"),
        target = c("Sample1",
                   "Sample1"),
        reference = c("Sample2", "Sample2"),
        n1 = 4:5,
        n2 = 4:5,
        statistic = c(W = 8, W = 12.5),
        p = c(1, 1),
        p_adj = c(1,
                  1),
        conf.low = c(-0.159036147906079,-8.12133146083903),
        conf.high = c(0.159036147906079, 8.12133146083903),
        method = c("Wilcoxon",
                   "Wilcoxon"),
        alternative = c("two.sided", "two.sided"),
        marker = c("ACTB",
                   "B2M")
      ),
      row.names = c(NA,-2L),
      class = c("tbl_df", "tbl",
                "data.frame")
    )

  expect_equal(dpa_markers[1:2, ], expected_result)

})

test_that("RunDPA fails with invalid input", {

  expect_error(
    dpa_markers <- RunDPA(polarization_table_merged),
    'argument "contrast_column" is missing, with no default'
  )
  expect_error(
    dpa_markers <- RunDPA(polarization_table_merged, reference = "Sample1"),
    'argument "contrast_column" is missing, with no default'
  )
  expect_error(
    dpa_markers <- RunDPA(polarization_table_merged, contrast_column = "sample", targets = "Sample1"),
    'argument "reference" is missing, with no default'
  )
  expect_error(
    dpa_markers <- RunDPA(polarization_table_merged, contrast_column = "Invalid", targets = "Sample1", reference = "Sample2")
  )
  expect_error(
    dpa_markers <- RunDPA(polarization_table_merged, contrast_column = "sample", targets = "Invalid", reference = "Sample2")
  )
  expect_error(
    dpa_markers <- RunDPA(polarization_table_merged, contrast_column = "sample", targets = "Sample1", reference = "Invalid")
  )
  expect_error(
    dpa_markers <- RunDPA(polarization_table_merged, contrast_column = "sample", targets = c("Sample1", "Sample2"), reference = "Sample1"),
    "targets is invalid\nall targets must be different from reference = 'Sample1'"
  )
  expect_error(
    dpa_markers <- RunDPA(polarization_table_merged, contrast_column = "sample", targets = "Sample2", reference = "Sample1", group_vars = "sample"),
    "contrast_column = 'sample' cannot be one of group_vars"
  )

  expect_error(
    dpa_markers <- RunDPA(polarization_table_merged, contrast_column = "sample", targets = "Sample2", reference = "Sample1", group_vars = "Invalid"),
    "group_vars is invalid\ngroup_vars must be a character vector with valid column names"
  )

})


if (TRUE) skip("Skipping parallel processing tests")

test_that("RunDPA can be parallelized", {

  # Sequential processing for reference
  expect_no_error(
    dpa_markers_seq <- RunDPA(polarization_table_merged, contrast_column = "sample", targets = "Sample2", reference = "Sample1", cl = 2)
  )

  # Using 2 threads. This will be ignored on Windows
  expect_no_error(
    dpa_markers_par <- RunDPA(polarization_table_merged, contrast_column = "sample", targets = "Sample2", reference = "Sample1", cl = 2)
  )

  expect_equal(dpa_markers_seq, dpa_markers_par)

  # Using a cluster object
  cl <- parallel::makeCluster(2)
  expect_no_error(
    dpa_markers <- RunDPA(polarization_table_merged, contrast_column = "sample", targets = "Sample2", reference = "Sample1", cl = cl)
  )
  parallel::stopCluster(cl)

})
