pxl_file <- system.file("extdata/PBMC_10_cells",
                        "Sample01_test.pxl",
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
seur_merged <- suppressWarnings(merge(seur1, seur2))

test_that("RunDPA works as expected on a data.frame", {

  expect_no_error(suppressWarnings(dpa_markers <- RunDPA(polarization_table_merged, contrast_column = "sample",
                        target = "Sample1", reference = "Sample2")))

  expect_true(all(dpa_markers$data_type == "morans_z"))
  expect_true(all(dpa_markers$target == "Sample1"))
  expect_true(all(dpa_markers$reference == "Sample2"))
  expect_true(all(dpa_markers$n1[1:4] == c(10, 9, 5, 7)))
  expect_true(all(dpa_markers$n2[1:4] == c(10, 9, 5, 7)))
  expect_true(all(dpa_markers$statistic[1:4] == c(50.0, 40.5, 12.5, 24.5)))
  expect_true(all(dpa_markers$method == "Wilcoxon"))
  expect_true(all(dpa_markers$marker[1:4] == c("ACTB", "CD11c", "CD14", "CD19")))
})

test_that("RunDPA works as expected on a Seurat object", {

  expect_no_error(suppressWarnings(dpa_markers <- RunDPA(seur_merged, contrast_column = "sample",
                                                         target = "Sample1", reference = "Sample2")))

  expect_true(all(dpa_markers$data_type == "morans_z"))
  expect_true(all(dpa_markers$target == "Sample1"))
  expect_true(all(dpa_markers$reference == "Sample2"))
  expect_true(all(dpa_markers$n1[1:4] == c(10, 9, 5, 7)))
  expect_true(all(dpa_markers$n2[1:4] == c(10, 9, 5, 7)))
  expect_true(all(dpa_markers$statistic[1:4] == c(50.0, 40.5, 12.5, 24.5)))
  expect_true(all(dpa_markers$method == "Wilcoxon"))
  expect_true(all(dpa_markers$marker[1:4] == c("ACTB", "CD11c", "CD14", "CD19")))
})

test_that("RunDPA fails with invalid input",  {
  expect_error(dpa_markers <- RunDPA(polarization_table_merged),
               'argument "contrast_column" is missing, with no default')
  expect_error(dpa_markers <- RunDPA(polarization_table_merged, contrast_column = "sample"),
               'argument "target" is missing, with no default')
  expect_error(dpa_markers <- RunDPA(polarization_table_merged, contrast_column = "sample", target = "Sample1"),
               'argument "reference" is missing, with no default')
  expect_error(dpa_markers <- RunDPA(polarization_table_merged, contrast_column = "Invalid", target = "Sample1", reference = "Sample2"),
               "'contrast_column' must be a valid column name")
  expect_error(dpa_markers <- RunDPA(polarization_table_merged, contrast_column = "sample", target = "Invalid", reference = "Sample2"),
               "'target' must be present in 'contrast_column' column")
  expect_error(dpa_markers <- RunDPA(polarization_table_merged, contrast_column = "sample", target = "Sample1", reference = "Invalid"),
               "'reference' must be present in 'contrast_column' column")
})
