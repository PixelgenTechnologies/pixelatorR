pxl_file <- system.file("extdata/five_cells",
                        "five_cells.pxl",
                        package = "pixelatorR")
seur_obj <- ReadMPX_Seurat(pxl_file)
pxl_file_new <- fs::path_temp("test.pxl")

# Load cell graphs and compute layouts
seur_obj <- seur_obj %>%
  LoadCellGraphs(cells = colnames(.)[1:2]) %>%
  ComputeLayout(dim = 3)

# Export layouts to a pxl file
seur_obj_subset <- subset(seur_obj, cells = colnames(seur_obj)[1:2])
WriteMPX_pxl_file(seur_obj_subset, pxl_file_new, export_layouts = TRUE, overwrite = TRUE)


test_that("ReadMPX_layouts works as expected", {

  # Note that the returned list is not ordered according with seur_obj_subset
  expect_no_error({layouts <- ReadMPX_layouts(filename = pxl_file_new)})
  expect_type(layouts, "list")
  expect_equal(names(layouts), "pmds_3d")
  expect_equal(dim(layouts[[1]][[1]]), c(3507, 4))
  expect_equal(layouts$pmds_3d[[1]] %>% head(),
               structure(
                 list(
                   name = c(
                     "CTAACGGCGAGAATCGTGAAGTGTA",
                     "TGTGGTATGTCGTTTAAGAATTAGT",
                     "TTTTGGTTAAGGCATAAATAGTGAG",
                     "CGGGCGTAGGATAGGTCTGGGTAGT",
                     "GGTAGGTCCACTCTTATAGGTTTGT",
                     "TATAAAGCAATAGTTGAAGGTATTT"
                   ),
                   x = c(
                     -36.9081349119502,
                     17.7141032144432,-17.1958361894217,
                     11.5713232744387,
                     -75.9335522318152,
                     82.3431503327363
                   ),
                   y = c(
                     -33.0424053966679,
                     -0.63606515583166,
                     -62.6655711057752,-3.26614274789496,
                     14.9023914761145,
                     7.63387064966907
                   ),
                   z = c(
                     10.865264588417,-34.610854988066,
                     -17.6882955492484,
                     68.6908927991003,
                     -24.7593144859087,
                     22.328929902093
                   )
                 ),
                 row.names = c(NA,-6L),
                 class = c("tbl_df",
                           "tbl", "data.frame")
               ))

})


test_that("ReadMPX_layouts fails with invalid input", {

  #Invalid file name
  expect_error({layouts <- ReadMPX_layouts(filename = "Invalid")}, "File 'Invalid' does not exist.")

  # Read from pxl file without pre-computed layouts
  expect_error({layouts <- ReadMPX_layouts(filename = pxl_file)})

})

