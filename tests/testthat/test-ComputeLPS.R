library(pixelatorR)
library(dplyr)
library(SeuratObject)

se <- ReadPNA_Seurat(minimal_pna_pxl_file(), verbose = FALSE) %>%
  LoadCellGraphs(cells = colnames(.)[1:2], verbose = FALSE)

cg_list <- CellGraphs(se)[1:2]
set.seed(123)
cg_small_list <- lapply(cg_list, function(cg) {
  nodes <- cg@cellgraph %>%
    pull(name) %>%
    sample(min(5000, length(cg@cellgraph)))
  subset(cg, nodes = nodes)
})
cgl <- CreateCellGraphList(cg_small_list)
cg_small <- cg_small_list[[1]]

test_that("CreateCellGraphList works as expected", {
  expect_s3_class(cgl, "CellGraphList")
  expect_type(cgl, "list")
  expect_equal(length(cgl), 2)
  expect_equal(names(cgl), names(cg_small_list))
  expect_s4_class(cgl[[1]], "CellGraph")
  expect_equal(length(cgl[1]), 1)
  expect_s3_class(cgl[1], "CellGraphList")
  expect_equal(length(lapply(cgl, identity)), 2)
  expect_s4_class(lapply(cgl, identity)[[1]], "CellGraph")
  mapped <- lapply.CellGraphList(cgl, identity)
  expect_s3_class(mapped, "CellGraphList")
  expect_equal(length(mapped), 2)
  expect_s4_class(mapped[[1]], "CellGraph")
})

test_that("CellGraphList subsetting, concatenation, and type checks work", {
  expect_s3_class(c(cgl[1], cgl[2]), "CellGraphList")
  expect_equal(length(c(cgl[1], cgl[2])), 2)
  expect_error(cgl[[1]] <- "Invalid")
  cgl_assigned <- cgl
  cgl_assigned[[1]] <- cg_small
  expect_s4_class(cgl_assigned[[1]], "CellGraph")
  expect_s3_class(cgl_assigned, "CellGraphList")
})

test_that("CellGraphList keeps NULL placeholders for unloaded graphs", {
  cgl_null <- CreateCellGraphList(list(a = cg_small, b = NULL))
  expect_s4_class(cgl_null[[1]], "CellGraph")
  expect_null(cgl_null[[2]])
  expect_equal(names(cgl_null), c("a", "b"))

  cgl_unloaded <- cgl
  cgl_unloaded[[1]] <- NULL
  expect_null(cgl_unloaded[[1]])
  expect_equal(length(cgl_unloaded), 2)
  expect_equal(names(cgl_unloaded), names(cgl))
  expect_s4_class(cgl_unloaded[[2]], "CellGraph")

  combined <- c(cgl[1], list(unloaded = NULL))
  expect_s3_class(combined, "CellGraphList")
  expect_equal(names(combined), c(names(cgl)[1], "unloaded"))
  expect_null(combined[["unloaded"]])
})

test_that("CreateCellGraphList fails when invalid input is provided", {
  expect_error(CreateCellGraphList("Invalid"))
  expect_error(CreateCellGraphList(list(cg_small)))
  expect_error(CreateCellGraphList(list(a = cg_small, a = cg_small)))
  expect_error(CreateCellGraphList(list(a = cg_small, b = "Invalid")))
  expect_s3_class(
    CreateCellGraphList(list(a = cg_small, b = NULL)),
    "CellGraphList"
  )
})

test_that("ComputeLPS.CellGraph stores a matrix as a layer", {
  expect_no_error(cg_lps <- ComputeLPS(cg_small, markers = c("B2M", "CD45")))
  expect_true("lps" %in% Layers(cg_lps))
  lps <- LayerData(cg_lps, layer = "lps")
  expect_equal(nrow(lps), length(Cells(cg_small)))
  expect_equal(ncol(lps), 2)
  expect_equal(colnames(lps), c("B2M", "CD45"))
  expect_equal(rownames(lps), Cells(cg_small))
})

test_that("ComputeLPS.CellGraph stores a vector in meta.data", {
  expect_no_error(cg_lps <- ComputeLPS(cg_small, markers = "B2M", mode = "all", name = "lps_b2m"))
  meta <- CellGraphData(cg_lps, slot = "meta.data")
  expect_true("lps_b2m" %in% colnames(meta))
  expect_equal(nrow(meta), length(Cells(cg_small)))
  expect_false("lps" %in% setdiff(Layers(cg_lps), "counts"))
})

test_that("ComputeLPS.CellGraphList works as expected", {
  expect_no_error(cgl_lps <- ComputeLPS(cgl, markers = "B2M", verbose = FALSE))
  expect_true("lps" %in% Layers(cgl_lps[[1]]))
  expect_true("lps" %in% Layers(cgl_lps[[2]]))
})

test_that("ComputeLPS.CellGraphList intersects missing markers", {
  expect_no_error(
    cgl_lps <- ComputeLPS(cgl, markers = c("B2M", "NotAMarker"), verbose = FALSE)
  )
  expect_equal(colnames(LayerData(cgl_lps[[1]], layer = "lps")), "B2M")
  expect_equal(colnames(LayerData(cgl_lps[[2]], layer = "lps")), "B2M")
})

test_that("ComputeLPS.CellGraphList warns when all markers are missing", {
  expect_warning(
    cgl_skip <- ComputeLPS(CreateCellGraphList(cgl[1]), markers = "NotAMarker", verbose = FALSE)
  )
  expect_false("lps" %in% setdiff(Layers(cgl_skip[[1]]), "counts"))
})

test_that("ComputeLPS.PNAAssay and Seurat only process loaded cellgraphs", {
  se_one <- se
  cgs <- lapply(CellGraphs(se_one), function(x) NULL)
  cgs[[1]] <- cg_small
  CellGraphs(se_one) <- cgs

  expect_no_error(pna_lps <- ComputeLPS(se_one[["PNA"]], markers = "B2M", verbose = FALSE))
  expect_true("lps" %in% Layers(CellGraphs(pna_lps)[[1]]))
  expect_null(CellGraphs(pna_lps)[[2]])

  expect_no_error(se_lps <- ComputeLPS(se_one, markers = "B2M", verbose = FALSE))
  expect_true("lps" %in% Layers(CellGraphs(se_lps)[[1]]))
  expect_null(CellGraphs(se_lps)[[2]])
})

test_that("ComputeLPS fails when invalid input is provided", {
  expect_error(ComputeLPS("Invalid"))
  expect_error(ComputeLPS(cg_small, markers = "Invalid"))
  expect_error(ComputeLPS(cg_small, mode = "Invalid"))
  expect_error(ComputeLPS(cg_small, method = "Invalid"))
})
