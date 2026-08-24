library(dplyr)
library(tidygraph)

se <- ReadPNA_Seurat(minimal_pna_pxl_file(), verbose = FALSE)
se <- LoadCellGraphs(se, cells = colnames(se)[4], add_layout = TRUE, verbose = FALSE)
cg <- CellGraphs(se)[[4]]
all_nodes <- cg@cellgraph %N>% pull(name)
seed1 <- all_nodes[1]

test_that("distance_from_node_set works as expected", {
  expect_no_error(out <- distance_from_node_set(cg, seed1))

  # Returns a CellGraph with the new node attribute
  expect_s4_class(out, "CellGraph")
  out_dist <- out@cellgraph %N>% as_tibble() %>% select(name, distance_from_seed)
  expect_s3_class(out_dist, "tbl_df")
  expect_equal(names(out_dist), c("name", "distance_from_seed"))
  expect_type(out_dist$distance_from_seed, "integer")
  expect_equal(nrow(out_dist), length(all_nodes))

  # The seed itself has distance 0 and is the only one
  expect_equal(sum(out_dist$distance_from_seed == 0L, na.rm = TRUE), 1L)
  expect_equal(
    out_dist$distance_from_seed[out_dist$name == seed1], 0L
  )

  # Data-level snippet of the first 10 rows
  expect_equal(
    out_dist %>% slice(1:10),
    structure(
      list(
        name = c(
          "55840301536286099-umi1", "3550955672846941-umi2",
          "20167079200998091-umi1", "33849376607686952-umi2",
          "56704973007062051-umi1", "70796622900050774-umi2",
          "57946909517865241-umi1", "50521531064131258-umi2",
          "15898272712678452-umi1", "55203010621091460-umi2"
        ),
        distance_from_seed = c(0L, 1L, 10L, 11L, 12L, 11L, 8L, 9L, 10L, 11L)
      ),
      row.names = c(NA, -10L),
      class = c("tbl_df", "tbl", "data.frame")
    )
  )

  # Counts of nodes at the first few distance levels
  dist_counts <- as.integer(table(out_dist$distance_from_seed)[1:10])
  expect_equal(
    dist_counts,
    c(1L, 8L, 37L, 188L, 404L, 1202L, 1808L, 4191L, 4551L, 9130L)
  )
})

test_that("distance_from_node_set handles multiple seeds", {
  seeds_multi <- all_nodes[1:5]
  out <- distance_from_node_set(cg, seeds_multi)
  d <- out@cellgraph %N>% pull(distance_from_seed)

  # Each seed is its own zero-distance node
  expect_equal(sum(d == 0L, na.rm = TRUE), length(seeds_multi))
})

test_that("distance_from_node_set respects max_iter", {
  out <- distance_from_node_set(cg, seed1, max_iter = 3L)
  d <- out@cellgraph %N>% pull(distance_from_seed)

  # No distance beyond max_iter
  expect_equal(max(d, na.rm = TRUE), 3L)

  # Counts at distances 0-3 match the uncapped run
  expect_equal(
    as.integer(table(d)),
    c(1L, 8L, 37L, 188L)
  )

  # All remaining nodes are unreached (NA)
  expect_equal(sum(is.na(d)), length(all_nodes) - (1L + 8L + 37L + 188L))
})

test_that("distance_from_node_set replaces an existing distance_from_seed column", {
  cg_pre <- cg
  cg_pre@cellgraph <- cg_pre@cellgraph %N>% mutate(distance_from_seed = 999L)
  out <- distance_from_node_set(cg_pre, seed1)
  d <- out@cellgraph %N>% pull(distance_from_seed)

  # The pre-existing column has been overwritten — no sentinel 999 values remain
  expect_false(any(d == 999L, na.rm = TRUE))
  expect_equal(sum(d == 0L, na.rm = TRUE), 1L)
})

test_that("distance_from_node_set is verbose when requested", {
  expect_message(
    distance_from_node_set(cg, seed1, max_iter = 2L, verbose = TRUE),
    regexp = "Iteration"
  )
})

test_that("distance_from_node_set fails with invalid input", {
  expect_error(distance_from_node_set("Invalid", seed1))
  expect_error(distance_from_node_set(cg, 1L))
  expect_error(distance_from_node_set(cg, character(0)))
  expect_error(
    distance_from_node_set(cg, "not_a_real_node"),
    regexp = "seed nodes must be present"
  )
})
