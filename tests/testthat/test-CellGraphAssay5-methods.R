options(Seurat.object.assay.version = "v5")

pxl_file <- system.file("extdata/five_cells", "five_cells.pxl", package = "pixelatorR")
seur_obj <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)
cg_assay5 <- seur_obj[["mpxCells"]]

# CellGraphs method
test_that("CellGraphs.CellGraphAssay5 getter/setter works as expected", {
  cg_list <- CellGraphs(cg_assay5)
  expect_equal(
    cg_list,
    list(
      RCVCMP0000217 = NULL,
      RCVCMP0000118 = NULL,
      RCVCMP0000487 = NULL,
      RCVCMP0000655 = NULL,
      RCVCMP0000263 = NULL
    )
  )
  CellGraphs(cg_assay5) <- cg_assay5@cellgraphs
  cg_list <- cg_assay5@cellgraphs
  expect_equal(
    cg_list,
    list(
      RCVCMP0000217 = NULL,
      RCVCMP0000118 = NULL,
      RCVCMP0000487 = NULL,
      RCVCMP0000655 = NULL,
      RCVCMP0000263 = NULL
    )
  )
})

test_that("CellGraphs.CellGraphAssay5 getter/setter fails when invalid input is provided", {
  expect_error(CellGraphs("Invalid input"), "no applicable method for 'CellGraphs'")
  expect_error(CellGraphs(cg_assay5) <- "Invalid input", "Invalid class 'character'")
  expect_error(CellGraphs(cg_assay5) <- setNames(cg_assay5@cellgraphs, nm = paste0(names(cg_assay5@cellgraphs), "invalid")))
  cgs <- cg_assay5@cellgraphs
  cgs[[1]] <- "Invalid"
  expect_error({CellGraphs(cg_assay5) <- cgs})
})

# RenameCells method
test_that("RenameCells.CellGraphAssay5 method works as expected", {
  cg_assay5_renamed <- RenameCells(cg_assay5,
                                  new.names = paste0("A_", colnames(cg_assay5)))
  expect_s4_class(cg_assay5_renamed, "CellGraphAssay5")
  expect_equal(colnames(cg_assay5_renamed), paste0("A_", colnames(cg_assay5)))
})

test_that("RenameCells.CellGraphAssay5 method fails when invalid input is provided", {
  expect_error(RenameCells(cg_assay5, new.names = "Invalid"), "'new.names' must be a character vector where")
})

# subset method
test_that("subset.CellGraphAssay5 works as expected", {
  expect_no_error(suppressWarnings(cg_assay5_subset <- subset(cg_assay5, cells = colnames(cg_assay5)[1:2])))
  expect_equal(ncol(cg_assay5_subset), 2)
  expect_equal(colnames(cg_assay5_subset), c("RCVCMP0000217", "RCVCMP0000118"))
})

# merge method
cg_assay5 <- seur_obj[["mpxCells"]]
test_that("merge.CellGraphAssay5 works as expected", {
  expect_no_error(cg_assay5_merged <- merge(cg_assay5, y = cg_assay5, add.cell.ids = c("A", "B")))
  expect_equal(ncol(cg_assay5_merged), 10)
  expect_equal(length(CellGraphs(cg_assay5_merged)), 10)
  expect_no_error(cg_assay5_merged <- merge(cg_assay5, y = list(cg_assay5, cg_assay5), add.cell.ids = c("A", "B", "C")))
  expect_equal(ncol(cg_assay5_merged), 15)
  expect_equal(length(CellGraphs(cg_assay5_merged)), 15)
  expect_no_error(cg_assay5_merged <- merge(cg_assay5, y = list(cg_assay5, cg_assay5), add.cell.ids = c("A", "B", "C")))
  expect_equal(colnames(cg_assay5_merged), c(paste0("A_", colnames(cg_assay5)),
                                            paste0("B_", colnames(cg_assay5)),
                                            paste0("C_", colnames(cg_assay5))))
  expect_no_error({cg_assay5_double_merged <- merge(cg_assay5_merged, cg_assay5_merged, add.cell.ids = c("A", "B"))})
})

test_that("merge.CellGraphAssay5 fails when invalid input is provided", {
  expect_error({cg_assay5_merged <- merge(cg_assay5, y = "Invalid")}, "'y' must be a 'CellGraphAssay5' object or a list of 'CellGraphAssay5' objects")
  expect_error({cg_assay5_merged <- merge(cg_assay5, y = list(cg_assay5, "Invalid"))}, "Element 2 in 'y' is not a 'CellGraphAssay5'")
  expect_no_error({cg_assay5_merged <- merge(cg_assay5, y = list(cg_assay5, cg_assay5), add.cell.ids = c("A", "B", "C"))})
  expect_no_error({cg_assay5_double_merged <- merge(cg_assay5_merged, cg_assay5_merged, add.cell.ids = c("A", "B"))})
})

# Show method
test_that("show.CellGraphAssay5 works as expected", {
  expect_no_error(capture.output(show(cg_assay5)))
})

# PolarizationScores getter/setter method
test_that("PolarizationScores.CellGraphAssay5 works as expected", {

  # Getter
  expect_no_error(pol <- PolarizationScores(cg_assay5))
  expect_s3_class(pol, "tbl_df")
  expect_equal(names(pol), c("morans_i", "morans_p_value", "morans_p_adjusted", "morans_z", "marker", "component"))
  expect_equal(dim(pol), c(380, 6))

  # Setter
  expect_no_error(PolarizationScores(cg_assay5) <- PolarizationScores(cg_assay5))
  expect_no_error(pol <- PolarizationScores(cg_assay5))
  expect_s3_class(pol, "tbl_df")
  expect_equal(names(pol), c("morans_i", "morans_p_value", "morans_p_adjusted", "morans_z", "marker", "component"))
  expect_equal(dim(pol), c(380, 6))
})

test_that("PolarizationScores.CellGraphAssay5 fails when invalid input is provided", {

  # Getter
  expect_error(pol <- PolarizationScores("invalid"), "no applicable method for")

  # Setter
  expect_error(PolarizationScores(cg_assay5) <- "invalid", "'polarization' must be a non-empty 'tbl_df' object")
})


# ColocalizationScores getter/setter method
test_that("ColocalizationScores.CellGraphAssay5 works as expected", {

  # Getter
  expect_no_error(coloc <- ColocalizationScores(cg_assay5))
  expect_s3_class(coloc, "tbl_df")
  expect_equal(names(coloc),
               c("marker_1", "marker_2", "pearson",
                 "pearson_mean", "pearson_stdev", "pearson_z",
                 "pearson_p_value", "pearson_p_value_adjusted",
                 "jaccard", "jaccard_mean", "jaccard_stdev",
                 "jaccard_z", "jaccard_p_value", "jaccard_p_value_adjusted",
                 "component"))
  expect_equal(dim(coloc), c(14649, 15))

  # Setter
  expect_no_error(ColocalizationScores(cg_assay5) <- ColocalizationScores(cg_assay5))
  expect_no_error(coloc <- ColocalizationScores(cg_assay5))
  expect_s3_class(coloc, "tbl_df")
  expect_equal(names(coloc),
               c("marker_1", "marker_2", "pearson",
                 "pearson_mean", "pearson_stdev", "pearson_z",
                 "pearson_p_value", "pearson_p_value_adjusted",
                 "jaccard", "jaccard_mean", "jaccard_stdev",
                 "jaccard_z", "jaccard_p_value", "jaccard_p_value_adjusted",
                 "component"))
  expect_equal(dim(coloc), c(14649, 15))
})

test_that("ColocalizationScores.CellGraphAssay5 fails when invalid input is provided", {

  # Getter
  expect_error(coloc <- ColocalizationScores("invalid"), "no applicable method for")

  # Setter
  expect_error(ColocalizationScores(cg_assay5) <- "invalid", "'colocalization' must be a non-empty 'tbl_df' object")
})

# Reset global option
options(Seurat.object.assay.version = "v3")
