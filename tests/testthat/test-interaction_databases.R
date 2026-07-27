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
    edges,
    structure(
      list(
        uniprot_a = c("A", "A"), uniprot_b = c("B", "C"), score = c(1, 3), evidence = c("x", "z"), resource = c(
          "test",
          "test"
        ),
        resource_version = c("v1", "v1")
      ),
      row.names = c(NA, -2L), class = c(
        "tbl_df",
        "tbl", "data.frame"
      )
    )
  )

  # Factor scores must decode via labels, not internal codes
  factor_edges <- normalise_interaction_edges(
    data.frame(
      a = "A",
      b = "B",
      score = factor("500", levels = c("100", "500")),
      stringsAsFactors = FALSE
    ),
    a_col = "a",
    b_col = "b",
    score_col = "score",
    resource = "test",
    resource_version = "v1"
  )
  expect_equal(factor_edges$score[[1]], 500)
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

  # version "latest" writes only one file; must not file.copy onto itself
  expect_no_error(
    latest_path <- save_interaction_database(
      edges = db$edges,
      database = "string",
      version = "latest",
      cache_dir = cache_dir
    )
  )
  expect_equal(basename(latest_path), "string_latest.rds")
  expect_true(file.exists(latest_path))
})

test_that("load_interaction_database fails with invalid input", {
  cache_dir <- tempfile("idb_missing_")
  dir.create(cache_dir)
  expect_error(load_interaction_database("string", "latest", cache_dir = cache_dir))
  expect_error(load_interaction_database("not_a_db", cache_dir = cache_dir))

  bad_path <- file.path(cache_dir, "string_latest.rds")
  saveRDS(list(not_edges = 1), bad_path)
  expect_error(load_interaction_database("string", "latest", cache_dir = cache_dir))

  saveRDS(
    list(edges = data.frame(x = 1), meta = list()),
    bad_path
  )
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
  expect_equal(
    out,
    structure(list(
      marker_1 = "A", marker_2 = "B", uniprot_a = "P12345",
      uniprot_b = "Q67890", in_db = TRUE, score = 500, evidence = "physical",
      resource = "string_physical", resource_version = "test"
    ), row.names = c(
      NA,
      -1L
    ), class = c("tbl_df", "tbl", "data.frame"))
  )

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

  # Overlapping UniProt mappings must not invent self-pairs via a Cartesian join
  overlap_map <- tibble(
    marker = c("A", "B", "A", "B"),
    uniprot_id = c("P12345", "Q67890", "Q67890", "P12345")
  )
  overlap_out <- extract_panel_interactions(
    markers = c("A", "B"),
    database = "string",
    marker_uniprot_map = overlap_map,
    score_threshold = 400,
    string_network = "physical",
    cache_dir = cache_dir
  )
  expect_equal(nrow(overlap_out), 1)
  expect_equal(overlap_out$marker_1[[1]], "A")
  expect_equal(overlap_out$marker_2[[1]], "B")
  expect_false(any(overlap_out$marker_1 == overlap_out$marker_2))
})

test_that("extract_panel_interactions keeps UniProt homodimers", {
  cache_dir <- tempfile("idb_homo_")
  dir.create(cache_dir)

  homo <- normalise_interaction_edges(
    data.frame(
      uniprot_a = c("P12345", "P12345"),
      uniprot_b = c("P12345", "Q67890"),
      score = c(900, 500),
      stringsAsFactors = FALSE
    ),
    score_col = "score",
    resource = "string_physical",
    resource_version = "test"
  )
  homo$evidence <- "physical"
  save_interaction_database(
    edges = homo,
    database = "string",
    version = "test",
    cache_dir = cache_dir
  )

  out <- extract_panel_interactions(
    markers = c("A", "B"),
    database = "string",
    marker_uniprot_map = marker_map,
    score_threshold = 400,
    string_network = "physical",
    cache_dir = cache_dir
  )
  expect_equal(nrow(out), 2)
  expect_true(any(out$marker_1 == out$marker_2 & out$uniprot_a == out$uniprot_b))
  self <- out %>% filter(marker_1 == marker_2)
  expect_equal(nrow(self), 1)
  expect_equal(self$marker_1[[1]], "A")
  expect_equal(self$uniprot_a[[1]], "P12345")
})

test_that("extract_panel_interactions does not invent hetero pairs from homodimers", {
  cache_dir <- tempfile("idb_homo_false_")
  dir.create(cache_dir)

  homo <- normalise_interaction_edges(
    data.frame(
      uniprot_a = "P12345",
      uniprot_b = "P12345",
      score = 900,
      stringsAsFactors = FALSE
    ),
    score_col = "score",
    resource = "string_physical",
    resource_version = "test"
  )
  homo$evidence <- "physical"
  save_interaction_database(
    edges = homo,
    database = "string",
    version = "test",
    cache_dir = cache_dir
  )

  shared_map <- tibble(
    marker = c("A", "B"),
    uniprot_id = c("P12345", "P12345")
  )
  out <- extract_panel_interactions(
    markers = c("A", "B"),
    database = "string",
    marker_uniprot_map = shared_map,
    score_threshold = 400,
    string_network = "physical",
    cache_dir = cache_dir
  )
  expect_equal(nrow(out), 2)
  expect_true(all(out$marker_1 == out$marker_2))
  expect_false(any(out$marker_1 != out$marker_2))
  expect_setequal(out$marker_1, c("A", "B"))
})

test_that("extract_panel_interactions keeps UniProt columns aligned with markers", {
  cache_dir <- tempfile("idb_align_")
  dir.create(cache_dir)

  edges <- normalise_interaction_edges(
    data.frame(
      uniprot_a = "P12345",
      uniprot_b = "Q67890",
      score = 700,
      stringsAsFactors = FALSE
    ),
    score_col = "score",
    resource = "string_physical",
    resource_version = "test"
  )
  edges$evidence <- "physical"
  save_interaction_database(
    edges = edges,
    database = "string",
    version = "test",
    cache_dir = cache_dir
  )

  # Join yields marker_1=B (P12345), marker_2=A (Q67890); lexical reorder to A,B
  # must swap UniProt columns with the markers.
  swapped_map <- tibble(
    marker = c("B", "A"),
    uniprot_id = c("P12345", "Q67890")
  )
  out <- extract_panel_interactions(
    markers = c("A", "B"),
    database = "string",
    marker_uniprot_map = swapped_map,
    score_threshold = 400,
    string_network = "physical",
    cache_dir = cache_dir
  )
  expect_equal(nrow(out), 1)
  expect_equal(out$marker_1[[1]], "A")
  expect_equal(out$marker_2[[1]], "B")
  expect_equal(out$uniprot_a[[1]], "Q67890")
  expect_equal(out$uniprot_b[[1]], "P12345")
})

test_that("create_marker_uniprot_map works as expected", {
  counts <- matrix(
    1:6,
    nrow = 3,
    dimnames = list(c("CD3e", "CD4", "CD8"), paste0("c", 1:2))
  )
  seur <- SeuratObject::CreateSeuratObject(counts = counts, assay = "PNA")
  seur[["PNA"]][[]] <- data.frame(
    uniprot_id = c("P07766", "P01730;P01732", "P01732"),
    row.names = c("CD3e", "CD4", "CD8"),
    stringsAsFactors = FALSE
  )

  expect_no_error(
    map <- create_marker_uniprot_map(seur, assay = "PNA")
  )
  expect_equal(names(map), c("marker", "uniprot_id"))
  expect_equal(nrow(map), 4)
  expect_true(all(c("CD3e", "CD4", "CD8") %in% map$marker))
  expect_equal(sum(map$marker == "CD4"), 2)
  expect_true(all(c("P01730", "P01732") %in% map$uniprot_id[map$marker == "CD4"]))

  expect_error(create_marker_uniprot_map(seur, assay = "missing"))
  seur2 <- SeuratObject::CreateSeuratObject(counts = counts, assay = "PNA")
  expect_error(create_marker_uniprot_map(seur2, assay = "PNA"))
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
  # Only the high-confidence pair; low-score and unscored pairs dropped
  expect_equal(nrow(db$edges), 1)
  expect_equal(db$edges$uniprot_a[[1]], "P12345")
  expect_equal(db$edges$uniprot_b[[1]], "Q67890")
  expect_true(all(db$edges$uniprot_a <= db$edges$uniprot_b))
  expect_equal(db$edges$score[[1]], 0.8)
})

test_that("build_alphafold_database score uses qualifying pDockQ2", {
  raw_dir <- tempfile("afdb_pdq_")
  cache_dir <- tempfile("afdb_pdq_cache_")
  dir.create(raw_dir)
  dir.create(cache_dir)

  # ipSAE below cutoff but present; kept only via pDockQ2
  utils::write.csv(
    data.frame(
      uniprot_ac_1 = "P12345",
      uniprot_ac_2 = "Q67890",
      ipSAE = 0.4,
      pDockQ2 = 0.5,
      stringsAsFactors = FALSE
    ),
    file.path(raw_dir, "afdb_api_complexes_panel.csv"),
    row.names = FALSE
  )

  build_alphafold_database(
    heterodimer_file = file.path(raw_dir, "missing.csv"),
    raw_dir = raw_dir,
    version = "test_af_pdq",
    cache_dir = cache_dir,
    ipsae_min = 0.6,
    pdockq2_min = 0.23
  )
  db <- load_interaction_database("alphafold", "latest", cache_dir = cache_dir)
  expect_equal(nrow(db$edges), 1)
  expect_equal(db$edges$score[[1]], 0.5)

  out <- extract_panel_interactions(
    markers = c("A", "B"),
    database = "alphafold",
    marker_uniprot_map = marker_map,
    score_threshold = 0.45,
    cache_dir = cache_dir
  )
  expect_equal(nrow(out), 1)
  out_hi <- extract_panel_interactions(
    markers = c("A", "B"),
    database = "alphafold",
    marker_uniprot_map = marker_map,
    score_threshold = 0.6,
    cache_dir = cache_dir
  )
  expect_equal(nrow(out_hi), 0)
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
  # 3 subunits -> 3 pairs; single-accession complex_2 -> 1 homomer edge
  expect_equal(nrow(db$edges), 4)
  expect_true(any(db$edges$uniprot_a == "A00001" & db$edges$uniprot_b == "A00001"))
})

test_that("build_corum_database fails with invalid input", {
  skip_if_not_installed("data.table")

  cache_dir <- tempfile("corum_bad_")
  dir.create(cache_dir)

  html_file <- tempfile(fileext = ".tsv")
  writeLines("<!DOCTYPE html><html></html>", html_file)
  expect_error(
    build_corum_database(corum_file = html_file, cache_dir = cache_dir)
  )
})

test_that("build_string_database prefers primary UniProt over secondary AC order", {
  skip_if_not_installed("data.table")

  raw_dir <- tempfile("string_raw_")
  cache_dir <- tempfile("string_cache_")
  dir.create(raw_dir)
  dir.create(cache_dir)

  # ITGAL-like: secondary AC O43746 precedes primary P20701; HGNC marks primary
  aliases_path <- file.path(raw_dir, "9606.protein.aliases.v12.0.txt.gz")
  phys_path <- file.path(raw_dir, "9606.protein.physical.links.v12.0.txt.gz")
  con <- gzfile(aliases_path, "wt")
  writeLines(
    c(
      "#string_protein_id\talias\tsource",
      "9606.ENSP00000349252\tO43746\tUniProt_AC",
      "9606.ENSP00000349252\tP20701\tUniProt_AC",
      "9606.ENSP00000349252\tP20701\tEnsembl_HGNC_uniprot_ids",
      "9606.ENSP00000380948\tP05107\tUniProt_AC",
      "9606.ENSP00000380948\tP05107\tEnsembl_HGNC_uniprot_ids"
    ),
    con
  )
  close(con)
  con <- gzfile(phys_path, "wt")
  writeLines(
    c(
      "protein1 protein2 combined_score",
      "9606.ENSP00000349252 9606.ENSP00000380948 999"
    ),
    con
  )
  close(con)

  path <- build_string_database(
    raw_dir = raw_dir,
    version = "12.0",
    cache_dir = cache_dir,
    species = 9606,
    include_full = FALSE
  )
  expect_true(file.exists(path))
  db <- load_interaction_database("string", "latest", cache_dir = cache_dir)
  pair_ids <- paste(db$edges$uniprot_a, db$edges$uniprot_b, sep = "-")
  expect_true("P05107-P20701" %in% pair_ids)
  expect_false("O43746-P05107" %in% pair_ids)

  out <- extract_panel_interactions(
    markers = c("CD11a", "CD18"),
    database = "string",
    marker_uniprot_map = tibble(
      marker = c("CD11a", "CD18"),
      uniprot_id = c("P20701", "P05107")
    ),
    score_threshold = 400,
    string_network = "physical",
    cache_dir = cache_dir
  )
  expect_equal(nrow(out), 1)
  expect_equal(out$marker_1[[1]], "CD11a")
  expect_equal(out$marker_2[[1]], "CD18")
  expect_equal(out$score[[1]], 999)
})

test_that(".download_if_missing reuses cached files and errors on bad URLs", {
  dest <- tempfile(fileext = ".txt")
  writeLines("cached", dest)
  expect_identical(
    .download_if_missing("https://example.invalid/should-not-fetch", dest),
    dest
  )
  expect_equal(readLines(dest, warn = FALSE), "cached")

  dl_dir <- tempfile("idb_dl_")
  dir.create(dl_dir)
  bad_dest <- file.path(dl_dir, "missing.bin")
  expect_error(
    .download_if_missing("https://example.invalid/missing.bin", bad_dest)
  )
  expect_false(file.exists(bad_dest))
  expect_equal(length(list.files(dl_dir, all.files = TRUE, no.. = TRUE)), 0)
})

test_that("build_alphafold_database accepts NVIDIA directional score columns", {
  raw_dir <- tempfile("afdb_nvda_")
  cache_dir <- tempfile("afdb_nvda_cache_")
  dir.create(raw_dir)
  dir.create(cache_dir)

  utils::write.csv(
    data.frame(
      uniprot_ac_1 = c("P12345", "P12345"),
      uniprot_ac_2 = c("Q67890", "P99999"),
      tax_id_1 = c(9606, 9606),
      tax_id_2 = c(9606, 10090),
      ipSAE_AB = c(0.9, 0.95),
      ipSAE_BA = c(0.7, 0.95),
      pDockQ2_AB = c(0.4, 0.5),
      pDockQ2_BA = c(0.3, 0.5),
      stringsAsFactors = FALSE
    ),
    file.path(raw_dir, "heterodimer_metadata.csv"),
    row.names = FALSE
  )

  path <- build_alphafold_database(
    heterodimer_file = file.path(raw_dir, "missing.csv"),
    raw_dir = raw_dir,
    version = "test_nvda",
    cache_dir = cache_dir,
    ipsae_min = 0.6,
    pdockq2_min = 0.23
  )
  db <- load_interaction_database("alphafold", "latest", cache_dir = cache_dir)
  expect_equal(nrow(db$edges), 1)
  expect_equal(db$edges$uniprot_a[[1]], "P12345")
  expect_equal(db$edges$uniprot_b[[1]], "Q67890")
  expect_true(file.exists(path))
})

test_that("build_alphafold_database parses homodimer NVIDIA schema", {
  raw_dir <- tempfile("afdb_homo_")
  cache_dir <- tempfile("afdb_homo_cache_")
  dir.create(raw_dir)
  dir.create(cache_dir)

  utils::write.csv(
    data.frame(
      uniprot_ac_1 = "P12345",
      uniprot_ac_2 = "Q67890",
      ipSAE = 0.8,
      pDockQ2 = 0.5,
      stringsAsFactors = FALSE
    ),
    file.path(raw_dir, "afdb_api_complexes_panel.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    data.frame(
      uniprotAccession = c("P99999", "A00001"),
      taxId = c(9606, 10090),
      ipSAE_AB = c(0.9, 0.95),
      ipSAE_BA = c(0.85, 0.95),
      pDockQ2_AB = c(0.4, 0.5),
      pDockQ2_BA = c(0.35, 0.5),
      stringsAsFactors = FALSE
    ),
    file.path(raw_dir, "afdb_homodimers_human_panel.csv"),
    row.names = FALSE
  )

  path <- build_alphafold_database(
    heterodimer_file = file.path(raw_dir, "missing.csv"),
    raw_dir = raw_dir,
    version = "test_af_homo",
    cache_dir = cache_dir,
    ipsae_min = 0.6,
    pdockq2_min = 0.23
  )
  expect_true(file.exists(path))
  db <- load_interaction_database("alphafold", "latest", cache_dir = cache_dir)
  # hetero P12345-Q67890 + human homodimer P99999-P99999; mouse A00001 dropped
  expect_equal(nrow(db$edges), 2)
  expect_true(any(db$edges$uniprot_a == "P99999" & db$edges$uniprot_b == "P99999"))
  expect_false(any(db$edges$uniprot_a == "A00001"))
})
