seur <- ReadMPX_Seurat(minimal_mpx_pxl_file())

seur_meta <- seur[[]]


test_that("sequencing_saturation works as expected", {
  expect_no_error(seqsat <-
    sequencing_saturation(
      graph_elements = seur_meta$edges,
      graph_reads = seur_meta$reads
    ))

  expect_equal(
    seqsat,
    c(
      80.2385969569762, 80.6342869229549, 81.2287328692286, 80.5187572104865,
      88.5155538178505
    )
  )

  expect_no_error(seqsat <-
    sequencing_saturation(
      graph_elements = seur_meta$vertices,
      graph_reads = seur_meta$reads
    ))

  expect_equal(
    seqsat,
    c(
      95.9017073454014, 95.3638094230871, 96.6128679798957, 95.5086128018459,
      98.5250840440276
    )
  )

  expect_error(sequencing_saturation(
    graph_elements = seur_meta$edges,
    graph_reads = NULL
  ))
  expect_error(suppressWarnings(
    sequencing_saturation(
      graph_elements = seur_meta$edges,
      graph_reads = "100"
    )
  ))
})


test_that("SequenceSaturationCurve works as expected", {
  edgelist <-
    ReadMPX_edgelist(system.file("extdata/five_cells", "five_cells.pxl", package = "pixelatorR")) %>%
    rename(umi1 = upia, umi2 = upib, read_count = count) %>%
    filter(component %in% c("RCVCMP0000118", "RCVCMP0000217"))

  set.seed(37)
  expect_no_error(seqsat <- SequenceSaturationCurve(edgelist,
    sample_fracs = c(1, 0.5),
    n_comps = 2L
  ))

  expect_equal(
    seqsat,
    structure(
      list(
        component = structure(c(1L, 1L, 2L, 2L),
          levels = c(
            "RCVCMP0000118",
            "RCVCMP0000217",
            "RCVCMP0000263",
            "RCVCMP0000487",
            "RCVCMP0000655"
          ), class = "factor"
        ),
        sample_size = c(75644, 37822, 60269, 30134), sample_frac = c(1, 0.5, 1, 0.5),
        graph_edges = c(
          14649L, 13416L,
          11910L, 10879L
        ),
        graph_proteins = c(3507L, 3430L, 2470L, 2418L),
        graph_reads = c(75644L, 37822L, 60269L, 30130L),
        graph_node_saturation = c(
          95.3638094230871,
          90.9312040611284,
          95.9017073454014,
          91.9747759707932
        ),
        graph_edge_saturation = c(
          80.6342869229549,
          64.5285812490085,
          80.2385969569762,
          63.8931297709924
        )
      ),
      class = c(
        "tbl_df",
        "tbl", "data.frame"
      ),
      row.names = c(NA, -4L)
    )
  )

  expect_error(SequenceSaturationCurve(edgelist,
    sample_fracs = c(1, 0.5),
    n_comps = 0
  ))
  expect_error(SequenceSaturationCurve(edgelist,
    sample_fracs = c(1, 0.5),
    n_comps = "2"
  ))
  expect_error(SequenceSaturationCurve(edgelist,
    sample_fracs = NULL
  ))
})


# approximate_edge_saturation
test_that("approximate_edge_saturation works as expected", {
  expect_no_error(db <- PixelDB$new(minimal_pna_pxl_file()))
  expect_no_error(edgesat <- approximate_edge_saturation(db))
  expect_equal(head(edgesat %>% collect() %>% arrange(component), 2),
               structure(
                 list(
                   component = c("0a45497c6bfbfb22", "2708240b908e2eba"),
                   edges = structure(c(
                     4.79312845656427e-319, 3.93463999035052e-319
                   ), class = "integer64"),
                   edge_saturation = c(0.751358719121778,
                                       0.742954902077026),
                   theoretical_max_edges = c(129118.086382752,
                                             107190.893790944)
                 ),
                 row.names = c(NA,-2L),
                 class = c("tbl_df",
                           "tbl", "data.frame")
               ))
  expect_error(approximate_edge_saturation("Invalid"))
  expect_error(approximate_edge_saturation(db, table_name = FALSE))
  expect_no_error(db$close())
})

# approximate_node_saturation
test_that("approximate_node_saturation works as expected", {
  expect_no_error(db <- PixelDB$new(minimal_pna_pxl_file()))
  expect_no_error(nodesat <- approximate_node_saturation(db))
  expect_equal(head(edgesat %>% collect() %>% arrange(component), 2),
               structure(
                 list(
                   component = c("0a45497c6bfbfb22", "2708240b908e2eba"),
                   nodes = c(43543L, 37665L),
                   node_saturation = c(0.885827322462419,
                                       0.886440632020902),
                   theoretical_max_nodes = c(49155.1783240997,
                                             42490.1551659828)
                 ),
                 row.names = c(NA,-2L),
                 class = c("tbl_df",
                           "tbl", "data.frame")
               ))
  expect_error(approximate_node_saturation("Invalid"))
  expect_error(approximate_node_saturation(db, table_name = FALSE))
  expect_no_error(db$close())
})


# approximate_saturation_curve
test_that("approximate_saturation_curve works as expected", {
  expect_no_error(db <- PixelDB$new(minimal_pna_pxl_file()))
  expect_no_error(sat_points <- approximate_saturation_curve(db, fracs = c(0.1, 0.2)))
  expect_equal(
    head(sat_points %>% collect() %>% arrange(component), 2),
    structure(
      list(
        component = c("0a45497c6bfbfb22", "0a45497c6bfbfb22"),
        p = c(0.1, 0.2),
        nodesat = c(0.546986224246965, 0.694512220887045),
        edgesat = c(0.183978283619467, 0.313676981034575),
        degree = c(1.7670057984165,
                   2.37274221920626),
        average_reads = c(29155, 58310)
      ),
      row.names = c(NA,-2L),
      class = c("tbl_df", "tbl", "data.frame")
    )
  )
  expect_error(approximate_saturation_curve("Invalid"))
  expect_error(approximate_saturation_curve(db, fracs = "Invalid"))
  expect_error(approximate_saturation_curve(db, detailed = "Invalid"))
  expect_error(approximate_saturation_curve(db, verbose = "Invalid"))
  expect_no_error(db$close())
})
