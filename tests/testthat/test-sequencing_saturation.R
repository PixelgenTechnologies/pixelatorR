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
