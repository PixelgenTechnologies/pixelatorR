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

  expect_equal(
    seqsat,
    structure(
      list(
        component = c(
          "2708240b908e2eba",
          "2708240b908e2eba",
          "c3c393e9a17c1981",
          "c3c393e9a17c1981"
        ),
        sample_size = c(
          229977,
          114988, 298516, 149258
        ),
        sample_frac = c(1, 0.5, 1, 0.5),
        graph_edges = c(
          79638L,
          57236L, 110657L, 78805L
        ),
        graph_proteins = c(
          37665L, 33064L,
          49351L, 43876L
        ),
        graph_reads = c(229977L, 113870L, 298516L, 148509L),
        graph_stability = c(1, 0.877844152396124, 1, 0.889059998784219),
        graph_node_saturation = c(
          83.6222752710054,
          70.9633792921753,
          83.4678878184084,
          70.4556626197739
        ),
        graph_edge_saturation = c(
          65.3713197406697,
          49.7356634758936,
          62.930965174396,
          46.9358759401787
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
  expect_equal(
    head(edgesat %>% collect() %>% arrange(component), 2),
    structure(
      list(
        component = c("0a45497c6bfbfb22", "2708240b908e2eba"),
        edges = c(97014L, 79638L),
        read_count = c(289856, 229977),
        edge_saturation = c(0.665302771031133, 0.653713197406697)
      ),
      row.names = c(NA, -2L),
      class = c("tbl_df", "tbl", "data.frame")
    )
  )
  expect_no_error(edgesat <- approximate_edge_saturation(db, components = "0a45497c6bfbfb22"))
  expect_equal(
    edgesat %>% collect(),
    structure(
      list(
        component = "0a45497c6bfbfb22",
        edges = 97014L,
        read_count = 289856,
        edge_saturation = 0.665302771031133
      ),
      class = c("tbl_df", "tbl", "data.frame"),
      row.names = c(NA, -1L)
    )
  )
  expect_no_error(edgesat_deg2 <- approximate_edge_saturation(
    db,
    components = "0a45497c6bfbfb22",
    min_degree = 2L
  ))
  expect_equal(
    edgesat_deg2 %>% collect(),
    structure(
      list(
        component = "0a45497c6bfbfb22",
        edges = 87202L,
        read_count = 270711,
        edge_saturation = 1 - (87202 / 270711)
      ),
      class = c("tbl_df", "tbl", "data.frame"),
      row.names = c(NA, -1L)
    )
  )
  expect_error(approximate_edge_saturation("Invalid"))
  expect_error(approximate_edge_saturation(db, table_name = FALSE))
  expect_error(approximate_edge_saturation(db, min_degree = 0L))
  expect_no_error(db$close())
})

# approximate_node_saturation
test_that("approximate_node_saturation works as expected", {
  expect_no_error(db <- PixelDB$new(minimal_pna_pxl_file()))
  expect_no_error(nodesat <- approximate_node_saturation(db))
  expect_equal(
    head(nodesat %>% collect() %>% arrange(component), 2),
    structure(
      list(
        component = c("0a45497c6bfbfb22", "2708240b908e2eba"),
        nodes = c(43543L, 37665L),
        read_count = c(579712, 459954),
        node_saturation = c(0.924888565356591, 0.918111376355027)
      ),
      row.names = c(NA, -2L),
      class = c("tbl_df", "tbl", "data.frame")
    )
  )
  expect_no_error(nodesat <- approximate_node_saturation(db, components = "0a45497c6bfbfb22"))
  expect_equal(
    nodesat %>% collect(),
    structure(
      list(
        component = "0a45497c6bfbfb22",
        nodes = 43543L,
        read_count = 579712,
        node_saturation = 0.924888565356591
      ),
      class = c("tbl_df", "tbl", "data.frame"),
      row.names = c(NA, -1L)
    )
  )
  expect_no_error(nodesat_m1 <- approximate_node_saturation(
    db,
    components = "0a45497c6bfbfb22",
    node_reads_multiplier = 1
  ))
  expect_equal(
    nodesat_m1 %>% collect(),
    structure(
      list(
        component = "0a45497c6bfbfb22",
        nodes = 43543L,
        read_count = 289856,
        node_saturation = 1 - (43543 / 289856)
      ),
      class = c("tbl_df", "tbl", "data.frame"),
      row.names = c(NA, -1L)
    )
  )
  expect_no_error(nodesat_deg2 <- approximate_node_saturation(
    db,
    components = "0a45497c6bfbfb22",
    min_degree = 2L
  ))
  expect_equal(
    nodesat_deg2 %>% collect(),
    structure(
      list(
        component = "0a45497c6bfbfb22",
        nodes = 33731L,
        read_count = 541422,
        node_saturation = 1 - (33731 / 541422)
      ),
      class = c("tbl_df", "tbl", "data.frame"),
      row.names = c(NA, -1L)
    )
  )
  expect_error(approximate_node_saturation("Invalid"))
  expect_error(approximate_node_saturation(db, table_name = FALSE))
  expect_error(approximate_node_saturation(db, node_reads_multiplier = 3))
  expect_error(approximate_node_saturation(db, min_degree = 0L))
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
        read_count = c(28986, 57971),
        p = c(0.1, 0.2),
        nodes = c(26393.2252179469, 33182.1442472842),
        edges = c(23754.9239169061, 40501.3715335435),
        degree = c(1.80007738506723, 2.44115456986227),
        node_saturation = c(0.544724604672136, 0.713803934318157),
        edge_saturation = c(0.180469056892771, 0.301351166384166)
      ),
      row.names = c(NA, -2L),
      class = c("tbl_df", "tbl", "data.frame")
    )
  )
  expect_no_error(sat_points_comp <- approximate_saturation_curve(
    db,
    fracs = c(0.1, 0.2),
    components = "0a45497c6bfbfb22",
    verbose = FALSE
  ))
  expect_equal(
    sat_points_comp,
    head(sat_points %>% collect() %>% arrange(component), 2)
  )
  expect_no_error(sat_points_deg2 <- approximate_saturation_curve(
    db,
    fracs = c(0.2),
    components = "0a45497c6bfbfb22",
    min_degree = 2L,
    verbose = FALSE
  ))
  expect_equal(sat_points_deg2$read_count, 57971)
  expect_lt(sat_points_deg2$nodes, sat_points_comp$nodes[2])
  expect_gt(sat_points_deg2$node_saturation, sat_points_comp$node_saturation[2])
  expect_error(approximate_saturation_curve("Invalid"))
  expect_error(approximate_saturation_curve(db, fracs = "Invalid"))
  expect_error(approximate_saturation_curve(db, detailed = "Invalid"))
  expect_error(approximate_saturation_curve(db, verbose = "Invalid"))
  expect_error(approximate_saturation_curve(db, min_degree = 0L))
  expect_error(approximate_saturation_curve(db, min_degree = 3L))
  expect_no_error(db$close())
})

# marginal_saturation and marginal_saturation_curve
test_that("marginal_saturation and marginal_saturation_curve work as expected", {
  expect_no_error(db <- PixelDB$new(minimal_pna_pxl_file()))
  expect_no_error(msat <- marginal_saturation(db, components = "0a45497c6bfbfb22") %>% collect())
  expect_s3_class(msat, "tbl_df")
  expect_equal(msat$read_count, 289856)
  expect_equal(msat$nodes, 43543L)
  expect_equal(msat$edges, 97014L)
  expect_equal(msat$degree, 2 * 97014 / 43543)
  expect_equal(msat$node_discovery_rate, 5694 / 289856)
  expect_equal(msat$edge_discovery_rate, 36446 / 289856)
  expect_equal(
    msat$degree_rate,
    (2 * (36446 / 289856) - (2 * 97014 / 43543) * (5694 / 289856)) / 43543
  )
  expect_equal(msat$marginal_node_saturation, 1 - (5694 / (2 * 289856)))
  expect_equal(msat$marginal_edge_saturation, 1 - (36446 / 289856))

  expect_no_error(msat_deg2 <- marginal_saturation(
    db,
    components = "0a45497c6bfbfb22",
    min_degree = 2L
  ) %>% collect())
  expect_equal(msat_deg2$nodes, 33731L)
  expect_equal(msat_deg2$degree, 2 * 97014 / 33731)
  expect_equal(msat_deg2$node_discovery_rate, 3494 / 289856)
  expect_equal(
    msat_deg2$degree_rate,
    (2 * (36446 / 289856) - (2 * 97014 / 33731) * (3494 / 289856)) / 33731
  )

  # Curve at p = 1.0 should match full-depth marginal_saturation()
  expect_no_error(mcurve_1 <- marginal_saturation_curve(
    db,
    fracs = c(0.2, 1.0),
    components = "0a45497c6bfbfb22",
    min_degree = 1L,
    verbose = FALSE
  ))
  expect_equal(mcurve_1$degree[2], msat$degree)
  expect_equal(mcurve_1$node_discovery_rate[2], msat$node_discovery_rate)
  expect_equal(mcurve_1$edge_discovery_rate[2], msat$edge_discovery_rate)
  expect_equal(mcurve_1$degree_rate[2], msat$degree_rate)
  expect_equal(mcurve_1$marginal_node_saturation[2], msat$marginal_node_saturation)
  expect_equal(mcurve_1$marginal_edge_saturation[2], msat$marginal_edge_saturation)

  expect_no_error(mcurve_2 <- marginal_saturation_curve(
    db,
    fracs = c(0.2, 1.0),
    components = "0a45497c6bfbfb22",
    min_degree = 2L,
    verbose = FALSE
  ))
  expect_equal(mcurve_2$nodes[2], as.numeric(msat_deg2$nodes))
  expect_equal(mcurve_2$degree[2], msat_deg2$degree)
  expect_equal(mcurve_2$node_discovery_rate[2], msat_deg2$node_discovery_rate)
  expect_equal(mcurve_2$degree_rate[2], msat_deg2$degree_rate)
  expect_equal(mcurve_2$marginal_node_saturation[2], msat_deg2$marginal_node_saturation)

  # Numerical derivative check for degree_rate at p = 0.5
  eps <- 1e-5
  mcurve_num <- marginal_saturation_curve(
    db,
    fracs = c(0.5 - eps, 0.5, 0.5 + eps),
    components = "0a45497c6bfbfb22",
    min_degree = 1L,
    verbose = FALSE
  )
  num_deg_rate <- (mcurve_num$degree[3] - mcurve_num$degree[1]) / (2 * eps * msat$read_count)
  expect_equal(mcurve_num$degree_rate[2], num_deg_rate, tolerance = 1e-6)

  expect_error(marginal_saturation("Invalid"))
  expect_error(marginal_saturation(db, min_degree = 3L))
  expect_error(marginal_saturation_curve("Invalid"))
  expect_error(marginal_saturation_curve(db, fracs = c(0, 0.5)))
  expect_error(marginal_saturation_curve(db, min_degree = 3L))

  # Uncollected lazy results should not be overwritten by later calls
  msat_lazy_1 <- marginal_saturation(db, components = "0a45497c6bfbfb22", min_degree = 1L)
  msat_lazy_2 <- marginal_saturation(db, components = "2708240b908e2eba", min_degree = 2L)
  expect_equal(msat_lazy_1 %>% collect(), msat)
  expect_equal((msat_lazy_2 %>% collect())$component, "2708240b908e2eba")
  expect_equal((msat_lazy_2 %>% collect())$nodes, 28944L)

  # Mixed use with components_edgelist(lazy = TRUE) on the same connection
  expect_no_error(lazy_el <- db$components_edgelist("0a45497c6bfbfb22", lazy = TRUE))
  expect_no_error(marginal_saturation(db, components = "0a45497c6bfbfb22") %>% collect())
  expect_no_error(db$components_edgelist("0a45497c6bfbfb22", lazy = TRUE))
  expect_no_error(approximate_edge_saturation(db, components = "0a45497c6bfbfb22") %>% collect())
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
if (!utils::compareVersion(as.character(duckdb_v), "1.5.0") > 0) {
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
    expect_equal(dim(lcc_pts), c(1, 4))
    expect_error(lcc_curve("Invalid"))
    expect_error(lcc_curve(minimal_pna_pxl_file(), "Invalid"))
    expect_error(lcc_curve(minimal_pna_pxl_file(), "0a45497c6bfbfb22", fracs = "Invalid"))
    expect_error(lcc_curve(minimal_pna_pxl_file(), "0a45497c6bfbfb22", outdir = FALSE))
  })
}
