library(dplyr)
library(tibble)

# Tiny STRING-like edge cache for extract / load tests
.setup_string_cache <- function(cache_dir) {
  phys <- normalise_interaction_edges(
    data.frame(
      uniprot_a = c("P12345", "P12345"),
      uniprot_b = c("Q67890", "P99999"),
      score = c(500, 200),
      stringsAsFactors = FALSE
    ),
    score_col = "score",
    resource = "string_physical",
    resource_version = "test"
  )
  phys$evidence <- "physical"

  full <- normalise_interaction_edges(
    data.frame(
      uniprot_a = "Q11111",
      uniprot_b = "R22222",
      score = 800,
      stringsAsFactors = FALSE
    ),
    score_col = "score",
    resource = "string_full",
    resource_version = "test"
  )
  full$evidence <- "full"

  save_interaction_database(
    edges = bind_rows(phys, full),
    database = "string",
    version = "test",
    cache_dir = cache_dir,
    source_url = "https://example.com",
    license = "CC BY 4.0"
  )
}

marker_map <- tibble(
  marker = c("A", "B", "C", "D"),
  uniprot_id = c("P12345", "Q67890", "P99999", "Q11111")
)

test_that("normalise_interaction_edges works as expected", {
  edges <- normalise_interaction_edges(
    data.frame(
      a = c("B", "A", "A", ""),
      b = c("A", "B", "C", "D"),
      score = c(1, 1, 3, 4),
      evid = c("x", "x", "z", "w"),
      stringsAsFactors = FALSE
    ),
    a_col = "a",
    b_col = "b",
    score_col = "score",
    evidence_col = "evid",
    resource = "test",
    resource_version = "v1"
  )

  expect_equal(
    edges %>%
      select(-built_at),
    structure(list(uniprot_a = c("A", "A"), uniprot_b = c("B", "C"
    ), score = c(1, 3), evidence = c("x", "z"), resource = c("test",
                                                             "test"),
    resource_version = c("v1", "v1")),
    row.names = c(NA, -2L), class = c("tbl_df",
                                      "tbl", "data.frame"))
  )
})

test_that("normalise_interaction_edges fails with invalid input", {
  expect_error(
    normalise_interaction_edges(
      data.frame(x = 1, y = 2),
      a_col = "a",
      b_col = "b",
      resource = "test",
      resource_version = "v1"
    )
  )
})

test_that("interaction_database_cache_dir works as expected", {
  expect_true(grepl("interaction_databases$", interaction_database_cache_dir()))
})

test_that("save_interaction_database and load_interaction_database work as expected", {
  cache_dir <- tempfile("idb_")
  dir.create(cache_dir)

  expect_no_error(path <- .setup_string_cache(cache_dir))
  expect_true(file.exists(path))
  expect_true(file.exists(file.path(cache_dir, "string_latest.rds")))

  expect_no_error(db <- load_interaction_database("string", "latest", cache_dir = cache_dir))
  expect_true(is.list(db))
  expect_true(inherits(db$edges, "tbl_df") || is.data.frame(db$edges))
  expect_equal(db$meta$database, "string")
  expect_equal(nrow(db$edges), 3)

  expect_no_error(db_ver <- load_interaction_database("string", "test", cache_dir = cache_dir))
  expect_equal(db_ver$meta$version, "test")
})

test_that("load_interaction_database fails with invalid input", {
  cache_dir <- tempfile("idb_missing_")
  dir.create(cache_dir)
  expect_error(load_interaction_database("string", "latest", cache_dir = cache_dir))
  expect_error(load_interaction_database("not_a_db", cache_dir = cache_dir))

  bad_path <- file.path(cache_dir, "string_latest.rds")
  saveRDS(list(not_edges = 1), bad_path)
  expect_error(load_interaction_database("string", "latest", cache_dir = cache_dir))
})

test_that("extract_panel_interactions works as expected", {
  cache_dir <- tempfile("idb_extract_")
  dir.create(cache_dir)
  .setup_string_cache(cache_dir)

  expect_no_error(
    out <- extract_panel_interactions(
      markers = c("A", "B", "C"),
      database = "string",
      marker_uniprot_map = marker_map,
      score_threshold = 400,
      string_network = "physical",
      cache_dir = cache_dir
    )
  )
  expect_equal(nrow(out), 1)
  expect_equal(out$marker_1, "A")
  expect_equal(out$marker_2, "B")
  expect_true(out$in_db)
  expect_equal(out$score, 500)

  # Lower threshold keeps the weaker physical edge A-C
  expect_no_error(
    out_low <- extract_panel_interactions(
      markers = c("A", "B", "C"),
      database = "string",
      marker_uniprot_map = marker_map,
      score_threshold = 100,
      string_network = "physical",
      cache_dir = cache_dir
    )
  )
  expect_equal(nrow(out_low), 2)

  # Markers with no UniProt overlap return an empty table
  expect_no_error(
    empty <- extract_panel_interactions(
      markers = "Z",
      database = "string",
      marker_uniprot_map = tibble(marker = "Z", uniprot_id = "P00000"),
      cache_dir = cache_dir
    )
  )
  expect_equal(nrow(empty), 0)
  expect_equal(
    names(empty),
    c(
      "marker_1", "marker_2", "uniprot_a", "uniprot_b", "in_db",
      "score", "evidence", "resource", "resource_version"
    )
  )
})

test_that("extract_panel_interactions fails with invalid input", {
  cache_dir <- tempfile("idb_extract_bad_")
  dir.create(cache_dir)
  .setup_string_cache(cache_dir)

  expect_error(
    extract_panel_interactions(
      markers = c("A", "B"),
      database = "string",
      cache_dir = cache_dir
    )
  )
  expect_error(
    extract_panel_interactions(
      markers = c("A", "B"),
      database = "string",
      marker_uniprot_map = tibble(x = 1, y = 2),
      cache_dir = cache_dir
    )
  )
  expect_error(
    extract_panel_interactions(
      markers = c("A", "B"),
      database = "not_a_db",
      marker_uniprot_map = marker_map,
      cache_dir = cache_dir
    )
  )
  expect_error(
    extract_panel_interactions(
      markers = c(1, 2),
      database = "string",
      marker_uniprot_map = marker_map,
      cache_dir = cache_dir
    )
  )
})

test_that("build_alphafold_database works as expected", {
  raw_dir <- tempfile("afdb_raw_")
  cache_dir <- tempfile("afdb_cache_")
  dir.create(raw_dir)
  dir.create(cache_dir)

  utils::write.csv(
    data.frame(
      uniprot_ac_1 = c("P12345", "P12345", "A00001"),
      uniprot_ac_2 = c("Q67890", "P99999", "A00002"),
      ipSAE = c(0.8, 0.1, NA_real_),
      pDockQ2 = c(0.5, 0.1, NA_real_),
      stringsAsFactors = FALSE
    ),
    file.path(raw_dir, "afdb_api_complexes_panel.csv"),
    row.names = FALSE
  )

  expect_no_error(
    path <- build_alphafold_database(
      heterodimer_file = file.path(raw_dir, "missing.csv"),
      raw_dir = raw_dir,
      version = "test_af",
      cache_dir = cache_dir,
      ipsae_min = 0.6,
      pdockq2_min = 0.23
    )
  )
  expect_true(file.exists(path))

  db <- load_interaction_database("alphafold", "latest", cache_dir = cache_dir)
  # High-confidence pair + NA-score prefiltered pair; low-score pair dropped
  expect_equal(nrow(db$edges), 2)
  expect_true(all(db$edges$uniprot_a <= db$edges$uniprot_b))
})

test_that("build_corum_database works as expected", {
  skip_if_not_installed("data.table")

  raw_dir <- tempfile("corum_raw_")
  cache_dir <- tempfile("corum_cache_")
  dir.create(raw_dir)
  dir.create(cache_dir)
  corum_file <- file.path(raw_dir, "omnipath_complexes_corum.tsv")

  utils::write.table(
    data.frame(
      name = c("complex_1", "complex_2"),
      components = c("P12345_Q67890_P99999", "A00001"),
      stringsAsFactors = FALSE
    ),
    corum_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  expect_no_error(
    path <- build_corum_database(
      corum_file = corum_file,
      version = "test_corum",
      cache_dir = cache_dir
    )
  )
  expect_true(file.exists(path))

  db <- load_interaction_database("corum", "latest", cache_dir = cache_dir)
  # 3 subunits -> 3 undirected pairs
  expect_equal(nrow(db$edges), 3)
})

test_that("build_corum_database fails with invalid input", {
  skip_if_not_installed("data.table")

  cache_dir <- tempfile("corum_bad_")
  dir.create(cache_dir)
  expect_error(
    build_corum_database(
      corum_file = file.path(tempdir(), "missing_corum.tsv"),
      cache_dir = cache_dir
    )
  )

  html_file <- tempfile(fileext = ".tsv")
  writeLines("<!DOCTYPE html><html></html>", html_file)
  expect_error(
    build_corum_database(corum_file = html_file, cache_dir = cache_dir)
  )
})
