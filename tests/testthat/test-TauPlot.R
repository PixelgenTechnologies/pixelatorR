pxl_file <- system.file("extdata/five_cells",
                        "five_cells.pxl",
                        package = "pixelatorR")
seur_obj <- ReadMPX_Seurat(pxl_file)

test_that("TauPlot works for Seurat objects", {
  expect_no_error({tau_plot <- TauPlot(seur_obj)})
  expect_s3_class(tau_plot, "ggplot")
  expect_equal(
    structure(list(x = ~tau, y = ~umi_per_upia, colour = ~tau_type), class = "uneval"),
    tau_plot$mapping
  )

  # With 0.18 PXL file format
  expect_no_error(
    seur_obj@meta.data <- seur_obj[[]] %>%
      rename(mean_molecules_per_a_pixel = mean_umi_per_upia) %>%
      select(-all_of("umi_per_upia"))
  )
  expect_no_error({tau_plot <- TauPlot(seur_obj)})
  expect_s3_class(tau_plot, "ggplot")
  expect_equal(
    structure(list(x = ~tau, y = ~mean_molecules_per_a_pixel, colour = ~tau_type), class = "uneval"),
    tau_plot$mapping
  )
})

test_that("TauPlot works for data.frame-like objects", {
  expect_no_error({tau_plot <- TauPlot(seur_obj[[]])})
  expect_s3_class(tau_plot, "ggplot")
})

test_that("TauPlot fails with invalid input", {
  expect_error(TauPlot("invalid"))
  expect_error(TauPlot(seur_obj, group_by = "invalid"))
})
