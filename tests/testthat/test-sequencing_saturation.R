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
    ReadPNA_edgelist(minimal_pna_pxl_file(), lazy = FALSE)

  set.seed(37)
  expect_no_error(seqsat <- SequenceSaturationCurve(edgelist,
    sample_fracs = c(1, 0.5),
    n_comps = 2L
  ))

  expect_equal(seqsat,
               structure(
                 list(
                   component = c(
                     "2708240b908e2eba",
                     "2708240b908e2eba",
                     "c3c393e9a17c1981",
                     "c3c393e9a17c1981"
                   ),
                   sample_size = c(229977,
                                   114988, 298516, 149258),
                   sample_frac = c(1, 0.5, 1, 0.5),
                   graph_edges = c(79638L,
                                   57236L, 110657L, 78805L),
                   graph_proteins = c(37665L, 33064L,
                                      49351L, 43876L),
                   graph_reads = c(229977L, 113870L, 298516L, 148509L),
                   graph_stability = c(1, 0.877844152396124, 1, 0.889059998784219),
                   graph_node_saturation = c(
                     91.8111376355027,
                     85.4816896460876,
                     91.7339439092042,
                     85.2278313098869
                   ),
                   graph_edge_saturation = c(
                     65.3713197406697,
                     49.7356634758936,
                     62.930965174396,
                     46.9358759401787
                   )
                 ),
                 class = c("tbl_df",
                           "tbl", "data.frame"),
                 row.names = c(NA,-4L)
               ))

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
                   edges = c(97014L, 79638L),
                   edge_saturation = c(0.751358719121778,
                                       0.742954902077026),
                   theoretical_max_edges = c(129118.086382752,
                                             107190.893790944)
                 ),
                 row.names = c(NA,-2L),
                 class = c("tbl_df",
                           "tbl", "data.frame")
               ))
  expect_no_error(edgesat <- approximate_edge_saturation(db, components = "0a45497c6bfbfb22"))
  expect_equal(edgesat %>% collect(),
               structure(
                 list(
                   component = "0a45497c6bfbfb22",
                   edges = 97014L,
                   edge_saturation = 0.751358719121778,
                   theoretical_max_edges = 129118.086382752
                 ),
                 class = c("tbl_df",
                           "tbl", "data.frame"),
                 row.names = c(NA,-1L)
               ))
  expect_error(approximate_edge_saturation("Invalid"))
  expect_error(approximate_edge_saturation(db, table_name = FALSE))
  expect_no_error(db$close())
})

# approximate_node_saturation
test_that("approximate_node_saturation works as expected", {
  expect_no_error(db <- PixelDB$new(minimal_pna_pxl_file()))
  expect_no_error(nodesat <- approximate_node_saturation(db))
  expect_equal(head(nodesat %>% collect() %>% arrange(component), 2),
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
  expect_no_error(nodesat <- approximate_node_saturation(db, components = "0a45497c6bfbfb22"))
  expect_equal(nodesat %>% collect(),
               structure(
                 list(
                   component = "0a45497c6bfbfb22",
                   nodes = 43543L,
                   node_saturation = 0.885827322462419,
                   theoretical_max_nodes = 49155.1783240997
                 ),
                 class = c("tbl_df",
                           "tbl", "data.frame"),
                 row.names = c(NA,-1L)
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
        degree = c(
          1.7670057984165,
          2.37274221920626
        ),
        average_reads = c(29155, 58310)
      ),
      row.names = c(NA, -2L),
      class = c("tbl_df", "tbl", "data.frame")
    )
  )
  expect_error(approximate_saturation_curve("Invalid"))
  expect_error(approximate_saturation_curve(db, fracs = "Invalid"))
  expect_error(approximate_saturation_curve(db, detailed = "Invalid"))
  expect_error(approximate_saturation_curve(db, verbose = "Invalid"))
  expect_no_error(db$close())
})

# downsample_to_parquet
test_that("downsample_to_parquet works as expected", {
  expect_no_error(pg_files <- downsample_to_parquet(minimal_pna_pxl_file(), components = "0a45497c6bfbfb22", fracs = 0.1, outdir = tempdir()))
  expect_s3_class(pg_files, "tbl_df")
  expect_equal(dim(pg_files), c(1, 2))
  expect_error(downsample_to_parquet("Invalid"))
  expect_error(downsample_to_parquet(minimal_pna_pxl_file(), "Invalid"))
  expect_error(downsample_to_parquet(minimal_pna_pxl_file(), components = "0a45497c6bfbfb22", fracs = "Invalid"))
  expect_error(downsample_to_parquet(minimal_pna_pxl_file(), components = "0a45497c6bfbfb22", outdir = FALSE))
})

duckdb_v <- packageVersion("duckdb")
if (!utils::compareVersion(as.character(duckdb_v), "1.3.2") > 0) {
  # lcc_sizes
  test_that("lcc_sizes works as expected", {
    expect_no_error(pg_files <- downsample_to_parquet(minimal_pna_pxl_file(), components = "0a45497c6bfbfb22", fracs = 0.1, outdir = tempdir()))
    expect_no_error(lcc <- lcc_sizes(pg_files))
    expect_s3_class(lcc, "tbl_df")
    expect_equal(dim(lcc), c(1, 3))
    expect_error(lcc_sizes("Invalid"))
    expect_error(lcc_sizes(pg_files, "Invalid"))
  })

  # lcc_curve
  test_that("lcc_curve works as expected", {
    set.seed(123)
    expect_no_error(lcc_pts <- lcc_curve(minimal_pna_pxl_file(), components = "0a45497c6bfbfb22", fracs = 0.1, outdir = tempdir()))
    expect_s3_class(lcc_pts, "tbl_df")
    expect_equal(dim(lcc_pts), c(1, 3))
    expect_error(lcc_curve("Invalid"))
    expect_error(lcc_curve(minimal_pna_pxl_file(), "Invalid"))
    expect_error(lcc_curve(minimal_pna_pxl_file(), "0a45497c6bfbfb22", fracs = "Invalid"))
    expect_error(lcc_curve(minimal_pna_pxl_file(), "0a45497c6bfbfb22", outdir = FALSE))
  })
}
