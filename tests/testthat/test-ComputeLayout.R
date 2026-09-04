library(dplyr)

for (assay_version in c("v3", "v5")) {
  options(Seurat.object.assay.version = assay_version)

  se <- ReadMPX_Seurat(minimal_mpx_pxl_file())

  se <- LoadCellGraphs(se, cells = colnames(se)[1:2])

  test_that("ComputeLayout works as expected", {
    # Seurat object
    expect_no_error(se <- se %>% ComputeLayout())

    # CellGraphAssay
    layout_method <- "pmds_3d"
    expect_no_error(cg_assay <- se[["mpxCells"]] %>% ComputeLayout(layout_method = "pmds"))
    expect_true(layout_method %in% names(cg_assay@cellgraphs[[1]]@layout))
    expect_equal(c("x", "y", "z"), colnames(cg_assay@cellgraphs[[1]]@layout[[layout_method]]))

    # Test with normalize_layout
    expect_no_error(cg_assay <- se[["mpxCells"]] %>%
      ComputeLayout(dim = 3, layout_method = "pmds", normalize_layout = TRUE))
    median_radius <- cg_assay@cellgraphs[[1]]@layout[["pmds_3d"]] %>%
      mutate(across(x:z, ~ .x^2)) %>%
      rowSums() %>%
      sqrt() %>%
      median()
    expect_equal(median_radius, 1)

    # Test with project_on_unit_sphere
    expect_no_error(cg_assay <- se[["mpxCells"]] %>%
      ComputeLayout(dim = 3, layout_method = "pmds", project_on_unit_sphere = TRUE))
    radii <- cg_assay@cellgraphs[[1]]@layout[[layout_method]] %>%
      mutate(across(x:z, ~ .x^2)) %>%
      rowSums() %>%
      sqrt() %>%
      max()
    expect_equal(unique(radii) %>% length(), 1)

    # Weighted pmds
    layout_method <- "wpmds"
    layout_method_3d <- "wpmds_3d"
    expect_no_error(se <- se %>% ComputeLayout(layout_method = layout_method, dim = 3))
    expect_true(layout_method_3d %in% names(CellGraphs(se)[[1]]@layout))
    expect_equal(c("x", "y", "z"), colnames(CellGraphs(se)[[1]]@layout[[layout_method_3d]]))
  })

  test_that("ComputeLayout works as expected with a custom layout function", {
    custom_layout_fkn <- graphlayouts::layout_with_pmds

    cg <- CellGraphs(se)[[colnames(se)[1]]]

    # tbl_graph
    expect_no_error(layout <- ComputeLayout(cg@cellgraph, custom_layout_function = custom_layout_fkn, custom_layout_function_args = list(pivots = 100, dim = 3)))
    expect_equal(dim(layout), c(2470, 3))

    # CellGraph
    expect_no_error(cg_layout <- ComputeLayout(cg, custom_layout_function = custom_layout_fkn, custom_layout_function_args = list(pivots = 100, dim = 3)))
    expect_true("custom" %in% names(cg_layout@layout))
    expect_equal(dim(cg_layout@layout[["custom"]]), c(2470, 3))

    # CellGraphAssay
    expect_no_error(cg_assay_layout <- ComputeLayout(se[["mpxCells"]], custom_layout_function = custom_layout_fkn, custom_layout_function_args = list(pivots = 100, dim = 3)))
    expect_true("custom" %in% names(CellGraphs(cg_assay_layout)[[1]]@layout))
    expect_equal(dim(CellGraphs(cg_assay_layout)[[1]]@layout[["custom"]]), c(2470, 3))

    # Seurat
    expect_no_error(se_layout <- ComputeLayout(se, custom_layout_function = custom_layout_fkn, custom_layout_function_args = list(pivots = 100, dim = 3)))
    expect_true("custom" %in% names(CellGraphs(se_layout)[[1]]@layout))
    expect_equal(dim(CellGraphs(se_layout)[[1]]@layout[["custom"]]), c(2470, 3))

    # Test with new layout name
    expect_no_error(se_layout <- ComputeLayout(se,
      custom_layout_function = custom_layout_fkn,
      custom_layout_function_args = list(pivots = 100, dim = 3), layout_name = "my_layout"
    ))
    expect_true("my_layout" %in% names(CellGraphs(se_layout)[[1]]@layout))
    expect_equal(dim(CellGraphs(se_layout)[[1]]@layout[["my_layout"]]), c(2470, 3))
  })

  test_that("custom layouts are stored like the built-in layouts", {
    custom_layout_fkn <- graphlayouts::layout_with_pmds

    cg <- CellGraphs(se)[[colnames(se)[1]]]
    node_names <- cg@cellgraph %>% pull(name)

    # Custom layout functions return a matrix, which is converted to the same
    # coordinate table that the built-in layout methods return
    layout <- ComputeLayout(cg@cellgraph,
      custom_layout_function = custom_layout_fkn,
      custom_layout_function_args = list(pivots = 100, dim = 3)
    )
    expect_s3_class(layout, "tbl_df")
    expect_equal(colnames(layout), c("x", "y", "z"))

    cg_layout <- ComputeLayout(cg,
      custom_layout_function = custom_layout_fkn,
      custom_layout_function_args = list(pivots = 100, dim = 3)
    )
    expect_s3_class(cg_layout@layout[["custom"]], "data.frame")
    expect_equal(colnames(cg_layout@layout[["custom"]]), c("x", "y", "z"))
    expect_equal(rownames(cg_layout@layout[["custom"]]), node_names)
  })

  test_that("custom layout coordinates are matched to nodes by name", {
    cg <- CellGraphs(se)[[colnames(se)[1]]]
    node_names <- cg@cellgraph %>% pull(name)
    coords <- graphlayouts::layout_with_pmds(cg@cellgraph, pivots = 100, dim = 3)

    # A custom function that returns its coordinates in a different order than
    # the graph nodes
    shuffled_fkn <- function(g, ...) coords[rev(seq_len(nrow(coords))), , drop = FALSE]
    cg_layout <- ComputeLayout(cg,
      custom_layout_function = shuffled_fkn,
      custom_layout_function_args = list(dim = 3)
    )
    expect_equal(rownames(cg_layout@layout[["custom"]]), node_names)
    expect_equal(
      as.numeric(as.matrix(cg_layout@layout[["custom"]])),
      as.numeric(coords[node_names, ])
    )

    # Coordinates without node names fall back to the graph node order
    unnamed_fkn <- function(g, ...) unname(coords)
    cg_layout <- ComputeLayout(cg,
      custom_layout_function = unnamed_fkn,
      custom_layout_function_args = list(dim = 3)
    )
    expect_equal(rownames(cg_layout@layout[["custom"]]), node_names)
    expect_equal(
      as.numeric(as.matrix(cg_layout@layout[["custom"]])),
      as.numeric(coords)
    )

    # Coordinates for nodes that are not in the graph are rejected
    unknown_nodes_fkn <- function(g, ...) {
      rownames(coords) <- paste0("unknown", seq_len(nrow(coords)))
      coords
    }
    expect_error(
      ComputeLayout(cg,
        custom_layout_function = unknown_nodes_fkn,
        custom_layout_function_args = list(dim = 3)
      ),
      "missing"
    )
  })


  test_that("ComputeLayout fails when invalid input is provided", {
    # Invalid object
    expect_error(ComputeLayout("Invalid"))
    e <- tryCatch(se[["mpxCells"]] %>% ComputeLayout(layout_method = "Invalid", verbose = FALSE), error = function(e) e, silent = TRUE)
    expect_s3_class(e, "simpleError")

    # Invalid combination of normalize_layout and project_on_unit_sphere
    expect_error(
      se[["mpxCells"]] %>% ComputeLayout(normalize_layout = TRUE, project_on_unit_sphere = TRUE, dim = 3)
    )

    # Invalid combination of dim and project_on_unit_sphere
    expect_error(
      se[["mpxCells"]] %>% ComputeLayout(dim = 2, project_on_unit_sphere = TRUE),
      "Projecting onto a unit sphere is only possible for 3D layouts"
    )

    custom_layout_fkn <- function(g) {
      return("Invalid")
    }

    cg <- CellGraphs(se)[[colnames(se)[1]]]

    # tbl_graph
    expect_error(layout <- ComputeLayout(cg@cellgraph, custom_layout_function = custom_layout_fkn))

    # CellGraph
    expect_error(cg_layout <- ComputeLayout(cg, custom_layout_function = custom_layout_fkn))

    # CellGraphAssay
    expect_error(cg_assay_layout <- ComputeLayout(se[["mpxCells"]], custom_layout_function = custom_layout_fkn))

    # Seurat
    expect_error(se_layout <- ComputeLayout(se, custom_layout_function = custom_layout_fkn))
  })
}

se <- ReadPNA_Seurat(minimal_pna_pxl_file()) %>%
  LoadCellGraphs(cells = colnames(.)[1])

test_that("ComputeLayout works with 'cpmds' option", {
  expect_no_error(se <- se %>% ComputeLayout(layout_method = "cpmds"))
  expect_true("cpmds_3d" %in% names(CellGraphs(se)[[1]]@layout))
  expect_equal(c("x", "y", "z"), colnames(CellGraphs(se)[[1]]@layout[["cpmds_3d"]]))
})
