pxl_file <- minimal_pna_pxl_file()

for (assay_version in c("v3", "v5")) {
  options(Seurat.object.assay.version = assay_version)

  expected_class <- ifelse(assay_version == "v3", "PNAAssay", "PNAAssay5")

  seur_obj <- suppressWarnings(ReadPNA_Seurat(pxl_file, verbose = FALSE))
  pna_assay <- seur_obj[["PNA"]]

  # CellGraphs method
  test_that("CellGraphs.PNAAssay getter/setter works as expected", {
    cg_list <- CellGraphs(pna_assay)
    expect_type(cg_list, "list")
    expect_equal(cg_list %>% length(), 5)
    CellGraphs(pna_assay) <- pna_assay@cellgraphs
    cg_list <- pna_assay@cellgraphs
    expect_equal(cg_list %>% length(), 5)
    CellGraphs(pna_assay) <- NULL
  })

  test_that("CellGraphs.PNAAssay getter/setter fails when invalid input is provided", {
    expect_error(CellGraphs("Invalid input"))
    expect_error(CellGraphs(pna_assay) <- "Invalid input")
    expect_error(CellGraphs(pna_assay) <-
      setNames(pna_assay@cellgraphs,
        nm = paste0(names(pna_assay@cellgraphs), "invalid")
      ))
    cgs <- pna_assay@cellgraphs
    cgs[[1]] <- "Invalid"
    expect_error(CellGraphs(pna_assay) <- cgs)
  })

  # RenameCells method
  test_that("RenameCells.PNAAssay method works as expected", {
    pna_assay_renamed <- RenameCells(pna_assay, new.names = paste0("A_", colnames(pna_assay)))
    expect_s4_class(pna_assay_renamed, expected_class)
    expect_equal(colnames(pna_assay_renamed), paste0("A_", colnames(pna_assay)))
  })

  test_that("RenameCells.PNAAssay method fails when invalid input is provided", {
    expect_error(RenameCells(pna_assay, new.names = "Invalid"))
  })

  # JoinLayers method
  if (assay_version == "v5") {
    test_that("JoinLayers.PNAAssay5 method works as expected", {
      expect_no_error(assay1 <- assay2 <- as(pna_assay, "Assay5"))
      colnames(assay1) <- paste0("A_", colnames(assay1))
      colnames(assay2) <- paste0("B_", colnames(assay2))
      expect_no_error(assay_merged <- merge(assay1, assay2))
      expect_no_error(pna_assay_merged <- as.PNAAssay5(assay_merged))
      expect_no_error(pna_assay_merged <- JoinLayers(pna_assay_merged))
      expect_equal(dim(LayerData(pna_assay_merged, layer = "counts")), c(158, 10))
    })
  }

  # subset method
  test_that("subset.PNAAssay works as expected", {
    if (assay_version == "v3") {
      pna_assay_subset <- subset(pna_assay, cells = colnames(pna_assay)[1:2])
    } else {
      expect_warning({
        pna_assay_subset <- subset(pna_assay, cells = colnames(pna_assay)[1:2])
      })
    }
    expect_equal(ncol(pna_assay_subset), 2)
    expect_equal(colnames(pna_assay_subset), c("0a45497c6bfbfb22", "2708240b908e2eba"))

    # Make sure that samples are correctly dropped
    expect_no_error(pna_assay_big <- merge(pna_assay, list(pna_assay, pna_assay), add.cell.ids = c("A", "B", "C")))
    if (assay_version == "v3") {
      expect_no_error(pna_assay_slice <- subset(pna_assay_big, cells = colnames(pna_assay_big)[1:2]))
    } else {
      expect_warning({
        pna_assay_slice <- subset(pna_assay_big, cells = colnames(pna_assay_big)[1:2])
      })
    }
    expect_true(nrow(FSMap(pna_assay_slice)) == 1)
  })

  # merge method
  pna_assay <- seur_obj[["PNA"]]
  test_that("merge.PNAAssay works as expected", {
    expect_no_error(pna_assay_merged <- merge(pna_assay, y = pna_assay, add.cell.ids = c("A", "B")))
    expect_equal(ncol(pna_assay_merged), 10)
    expect_equal(length(CellGraphs(pna_assay_merged)), 10)
    expect_no_error(pna_assay_merged <- merge(pna_assay, y = list(pna_assay, pna_assay), add.cell.ids = c("A", "B", "C")))
    expect_equal(ncol(pna_assay_merged), 15)
    expect_equal(length(CellGraphs(pna_assay_merged)), 15)
    expect_no_error(pna_assay_merged <-
      merge(pna_assay, y = list(pna_assay, pna_assay), add.cell.ids = c("A", "B", "C")))
    expect_equal(colnames(pna_assay_merged), c(
      paste0("A_", colnames(pna_assay)),
      paste0("B_", colnames(pna_assay)),
      paste0("C_", colnames(pna_assay))
    ))
    expect_no_error({
      pna_assay_merged <-
        merge(pna_assay, y = list(pna_assay, pna_assay), add.cell.ids = c("A", "B", "C"))
    })
    expect_no_error({
      pna_assay_double_merged <-
        merge(pna_assay_merged, pna_assay_merged, add.cell.ids = c("A", "B"))
    })
  })

  test_that("merge.PNAAssay fails when invalid input is provided", {
    expect_error(
      pna_assay_merged <- merge(pna_assay, y = "Invalid")
    )
    expect_error(
      pna_assay_merged <- merge(pna_assay, y = list(pna_assay, "Invalid"))
    )
    expect_no_error(
      pna_assay_merged <- merge(pna_assay, y = list(pna_assay, pna_assay), add.cell.ids = c("A", "B", "C"))
    )
    expect_no_error(
      pna_assay_double_merged <- merge(pna_assay_merged, pna_assay_merged, add.cell.ids = c("A", "B"))
    )
  })

  # Show method
  test_that("show.PNAAssay works without errors", {
    expect_no_error(capture.output(show(pna_assay)))
  })
}
