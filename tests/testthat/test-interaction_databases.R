library(dplyr)
library(tibble)

# Tiny STRING-like edge cache for extract / load tests
.setup_string_cache <- function(cache_dir) {
  phys <- normalise_interaction_edges(
    data.frame(
      uniprot_a = c("P12345", "P12345"),
      uniprot_b = c("Q67890", "P99999"),
      combined_score = c(500, 200),
      network = "physical",
      stringsAsFactors = FALSE
    ),
    score_cols = "combined_score",
    additional_cols = "network"
  )
  full <- normalise_interaction_edges(
    data.frame(
      uniprot_a = "Q11111",
      uniprot_b = "R22222",
      combined_score = 800,
      network = "full",
      stringsAsFactors = FALSE
    ),
    score_cols = "combined_score",
    additional_cols = "network"
  )
  save_interaction_database(
    edges = bind_rows(phys, full),
    database = "string",
    version = "test",
    cache_dir = cache_dir,
    score_columns = "combined_score",
    additional_columns = "network",
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
      primary_score = c(1, 1, 3, 4),
      evid = c("x", "x", "z", "w"),
      stringsAsFactors = FALSE
    ),
    a_col = "a",
    b_col = "b",
    score_cols = "primary_score",
    additional_cols = "evid"
  )
  expect_equal(
    edges,
    tibble(
      uniprot_a = c("A", "A"),
      uniprot_b = c("B", "C"),
      primary_score = c(1, 3),
      evid = c("x", "z")
    )
  )

  factor_edges <- normalise_interaction_edges(
    data.frame(
      a = "A", b = "B",
      score = factor("500", levels = c("100", "500")),
      stringsAsFactors = FALSE
    ),
    a_col = "a", b_col = "b", score_cols = "score"
  )
  expect_equal(factor_edges$score[[1]], 500)

  expect_error(
    normalise_interaction_edges(data.frame(x = 1, y = 2), a_col = "a", b_col = "b")
  )
})

test_that("save and load interaction database round-trip", {
  cache_dir <- tempfile("idb_")
  dir.create(cache_dir)
  .setup_string_cache(cache_dir)

  db <- load_interaction_database("string", "latest", cache_dir = cache_dir)
  expect_equal(
    db$edges %>% arrange(uniprot_a, uniprot_b, network),
    tibble(
      uniprot_a = c("P12345", "P12345", "Q11111"),
      uniprot_b = c("P99999", "Q67890", "R22222"),
      combined_score = c(200, 500, 800),
      network = c("physical", "physical", "full")
    )
  )
  expect_equal(db$meta$score_columns, "combined_score")
  expect_equal(db$meta$additional_columns, "network")
  expect_equal(db$meta$resource, "string")
  expect_false(
    any(c("evidence", "resource", "resource_version", "in_db") %in% names(db$edges))
  )

  expect_error(
    load_interaction_database("string", "latest", cache_dir = tempfile("missing_")),
    "build_if_missing = TRUE"
  )
})

test_that("build_string_database from gz fixture", {
  skip_if_not_installed("data.table")

  raw_dir <- tempfile("string_raw_")
  cache_dir <- tempfile("string_cache_")
  dir.create(raw_dir)
  dir.create(cache_dir)

  # Secondary AC O43746 precedes primary P20701; HGNC marks primary
  aliases_path <- file.path(raw_dir, "9606.protein.aliases.v12.0.txt.gz")
  phys_path <- file.path(raw_dir, "9606.protein.physical.links.v12.0.txt.gz")
  con <- gzfile(aliases_path, "wt")
  writeLines(
    c(
      "#string_protein_id\talias\tsource",
      "9606.ENSP1\tO43746\tUniProt_AC",
      "9606.ENSP1\tP20701\tUniProt_AC",
      "9606.ENSP1\tP20701\tEnsembl_HGNC_uniprot_ids",
      "9606.ENSP2\tP05107\tUniProt_AC",
      "9606.ENSP2\tP05107\tEnsembl_HGNC_uniprot_ids"
    ),
    con
  )
  close(con)
  con <- gzfile(phys_path, "wt")
  writeLines(
    c("protein1 protein2 combined_score", "9606.ENSP1 9606.ENSP2 999"),
    con
  )
  close(con)

  expect_message(
    db <- load_interaction_database(
      "string",
      version = "12.0",
      cache_dir = cache_dir,
      build_if_missing = TRUE,
      raw_dir = raw_dir,
      include_full = FALSE
    ),
    "Building interaction database"
  )
  expect_equal(
    head(db$edges %>% arrange(uniprot_a, uniprot_b), 1),
    tibble(
      uniprot_a = "P05107",
      uniprot_b = "P20701",
      combined_score = 999,
      network = "physical"
    )
  )
  expect_true(file.exists(aliases_path))

  expect_no_message(
    db2 <- load_interaction_database(
      "string",
      version = "12.0",
      cache_dir = cache_dir,
      build_if_missing = TRUE,
      raw_dir = raw_dir
    )
  )
  expect_equal(nrow(db2$edges), 1)
})

test_that("build_biogrid_database from tab3 fixture", {
  raw_dir <- tempfile("biogrid_raw_")
  cache_dir <- tempfile("biogrid_cache_")
  dir.create(raw_dir)
  dir.create(cache_dir)
  raw_file <- file.path(
    raw_dir,
    "BIOGRID-ORGANISM-Homo_sapiens-5.0.259.tab3.txt"
  )
  utils::write.table(
    data.frame(
      `SWISS-PROT Accessions Interactor A` = c(
        "P12345", "P11111", "P99999"
      ),
      `SWISS-PROT Accessions Interactor B` = c(
        "Q67890", "Q22222", "P99999"
      ),
      `Organism ID Interactor A` = c(9606, 9606, 9606),
      `Organism ID Interactor B` = c(9606, 9606, 10090),
      `Experimental System` = c(
        "Two-hybrid", "Synthetic Lethality", "Affinity Capture"
      ),
      `Experimental System Type` = c("physical", "genetic", "physical"),
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    raw_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  build_biogrid_database(
    raw_file = raw_file,
    version = "5.0.259",
    cache_dir = cache_dir
  )
  db <- load_interaction_database("biogrid", "5.0.259", cache_dir = cache_dir)
  expect_equal(
    db$edges,
    tibble(
      uniprot_a = "P12345",
      uniprot_b = "Q67890",
      `Experimental System` = "Two-hybrid"
    )
  )
  expect_equal(db$meta$additional_columns, "Experimental System")
  expect_error(
    build_biogrid_database(
      raw_file = raw_file,
      version = "4.4.238",
      cache_dir = cache_dir
    )
  )
})

test_that("build_alphafold_database merges panel and homodimer sources", {
  raw_dir <- tempfile("afdb_raw_")
  cache_dir <- tempfile("afdb_cache_")
  dir.create(raw_dir)
  dir.create(cache_dir)

  utils::write.csv(
    data.frame(
      uniprot_ac_1 = c("P12345", "P12345"),
      uniprot_ac_2 = c("Q67890", "P99999"),
      ipSAE = c(0.8, 0.1),
      pDockQ2 = c(0.5, 0.1),
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

  build_alphafold_database(
    heterodimer_file = file.path(raw_dir, "missing.csv"),
    raw_dir = raw_dir,
    version = "test_af",
    cache_dir = cache_dir
  )
  db <- load_interaction_database("alphafold", "latest", cache_dir = cache_dir)
  expect_equal(
    db$edges %>% arrange(uniprot_a, uniprot_b),
    tibble(
      uniprot_a = c("P12345", "P12345", "P99999"),
      uniprot_b = c("P99999", "Q67890", "P99999"),
      ipSAE = c(0.1, 0.8, 0.9),
      pDockQ2 = c(0.1, 0.5, 0.4)
    )
  )
  expect_equal(db$meta$score_columns, c("ipSAE", "pDockQ2"))
})

test_that("extract_panel_interactions on STRING cache", {
  cache_dir <- tempfile("idb_extract_")
  dir.create(cache_dir)
  .setup_string_cache(cache_dir)

  expect_equal(
    extract_panel_interactions(
      markers = c("A", "B", "C"),
      database = "string",
      marker_uniprot_map = marker_map,
      score_min = c(combined_score = 400),
      string_network = "physical",
      cache_dir = cache_dir
    ),
    tibble(
      marker_1 = "A", marker_2 = "B",
      uniprot_a = "P12345", uniprot_b = "Q67890",
      combined_score = 500, network = "physical"
    )
  )

  expect_equal(
    extract_panel_interactions(
      markers = c("A", "B", "C"),
      database = "string",
      marker_uniprot_map = marker_map,
      score_min = c(combined_score = 100),
      score_max = c(combined_score = 400),
      string_network = "physical",
      cache_dir = cache_dir
    ),
    tibble(
      marker_1 = "A", marker_2 = "C",
      uniprot_a = "P12345", uniprot_b = "P99999",
      combined_score = 200, network = "physical"
    )
  )

  empty <- extract_panel_interactions(
    markers = "Z",
    database = "string",
    marker_uniprot_map = tibble(marker = "Z", uniprot_id = "P00000"),
    cache_dir = cache_dir
  )
  expect_equal(nrow(empty), 0)
  expect_equal(
    names(empty),
    c("marker_1", "marker_2", "uniprot_a", "uniprot_b", "combined_score", "network")
  )

  # Max-score UniProt edge wins when multiple map to the same marker pair
  align_dir <- tempfile("idb_align_")
  dir.create(align_dir)
  save_interaction_database(
    edges = normalise_interaction_edges(
      data.frame(
        uniprot_a = c("P11111", "P22222"),
        uniprot_b = c("Q11111", "Q22222"),
        combined_score = c(100, 900),
        network = "physical",
        stringsAsFactors = FALSE
      ),
      score_cols = "combined_score",
      additional_cols = "network"
    ),
    database = "string",
    version = "test",
    cache_dir = align_dir,
    score_columns = "combined_score",
    additional_columns = "network"
  )
  expect_equal(
    extract_panel_interactions(
      markers = c("A", "B"),
      database = "string",
      marker_uniprot_map = tibble(
        marker = c("A", "A", "B", "B"),
        uniprot_id = c("P11111", "P22222", "Q11111", "Q22222")
      ),
      string_network = "physical",
      cache_dir = align_dir
    ),
    tibble(
      marker_1 = "A", marker_2 = "B",
      uniprot_a = "P22222", uniprot_b = "Q22222",
      combined_score = 900, network = "physical"
    )
  )

  # Homodimer self-pairs kept; shared accession does not invent hetero pairs
  homo_dir <- tempfile("idb_homo_")
  dir.create(homo_dir)
  save_interaction_database(
    edges = normalise_interaction_edges(
      data.frame(
        uniprot_a = "P12345",
        uniprot_b = "P12345",
        combined_score = 900,
        network = "physical",
        stringsAsFactors = FALSE
      ),
      score_cols = "combined_score",
      additional_cols = "network"
    ),
    database = "string",
    version = "test",
    cache_dir = homo_dir,
    score_columns = "combined_score",
    additional_columns = "network"
  )
  expect_equal(
    extract_panel_interactions(
      markers = c("A", "B"),
      database = "string",
      marker_uniprot_map = tibble(
        marker = c("A", "B"),
        uniprot_id = c("P12345", "P12345")
      ),
      score_min = c(combined_score = 400),
      string_network = "physical",
      cache_dir = homo_dir
    ) %>% arrange(marker_1),
    tibble(
      marker_1 = c("A", "B"),
      marker_2 = c("A", "B"),
      uniprot_a = c("P12345", "P12345"),
      uniprot_b = c("P12345", "P12345"),
      combined_score = c(900, 900),
      network = c("physical", "physical")
    )
  )

  # Lexical marker reorder swaps UniProt columns with markers
  swap_dir <- tempfile("idb_swap_")
  dir.create(swap_dir)
  save_interaction_database(
    edges = normalise_interaction_edges(
      data.frame(
        uniprot_a = "P12345",
        uniprot_b = "Q67890",
        combined_score = 700,
        network = "physical",
        stringsAsFactors = FALSE
      ),
      score_cols = "combined_score",
      additional_cols = "network"
    ),
    database = "string",
    version = "test",
    cache_dir = swap_dir,
    score_columns = "combined_score",
    additional_columns = "network"
  )
  expect_equal(
    extract_panel_interactions(
      markers = c("A", "B"),
      database = "string",
      marker_uniprot_map = tibble(
        marker = c("B", "A"),
        uniprot_id = c("P12345", "Q67890")
      ),
      score_min = c(combined_score = 400),
      string_network = "physical",
      cache_dir = swap_dir
    ),
    tibble(
      marker_1 = "A", marker_2 = "B",
      uniprot_a = "Q67890", uniprot_b = "P12345",
      combined_score = 700, network = "physical"
    )
  )
})

test_that("extract_panel_interactions filters multi-score AlphaFold edges", {
  raw_dir <- tempfile("afdb_combine_")
  cache_dir <- tempfile("afdb_combine_cache_")
  dir.create(raw_dir)
  dir.create(cache_dir)

  utils::write.csv(
    data.frame(
      uniprot_ac_1 = c("P12345", "P12345", "P99999"),
      uniprot_ac_2 = c("Q67890", "P99999", "Q67890"),
      ipSAE = c(0.8, 0.4, 0.7),
      pDockQ2 = c(0.1, 0.5, 0.4),
      stringsAsFactors = FALSE
    ),
    file.path(raw_dir, "afdb_api_complexes_panel.csv"),
    row.names = FALSE
  )
  build_alphafold_database(
    heterodimer_file = file.path(raw_dir, "missing.csv"),
    raw_dir = raw_dir,
    version = "test_combine",
    cache_dir = cache_dir
  )
  af_map <- tibble(
    marker = c("A", "B", "C"),
    uniprot_id = c("P12345", "Q67890", "P99999")
  )

  expect_equal(
    extract_panel_interactions(
      markers = c("A", "B", "C"),
      database = "alphafold",
      marker_uniprot_map = af_map,
      score_min = c(ipSAE = 0.6, pDockQ2 = 0.45),
      score_combine = "any",
      cache_dir = cache_dir
    ) %>%
      arrange(marker_1, marker_2) %>%
      select(marker_1, marker_2, ipSAE, pDockQ2),
    tibble(
      marker_1 = c("A", "A", "B"),
      marker_2 = c("B", "C", "C"),
      ipSAE = c(0.8, 0.4, 0.7),
      pDockQ2 = c(0.1, 0.5, 0.4)
    )
  )
  expect_equal(
    nrow(extract_panel_interactions(
      markers = c("A", "B", "C"),
      database = "alphafold",
      marker_uniprot_map = af_map,
      score_min = c(ipSAE = 0.6, pDockQ2 = 0.45),
      score_combine = "all",
      cache_dir = cache_dir
    )),
    0
  )
  expect_equal(
    extract_panel_interactions(
      markers = c("A", "B", "C"),
      database = "alphafold",
      marker_uniprot_map = af_map,
      score_min = c(ipSAE = 0.3, pDockQ2 = 0.45),
      score_max = c(ipSAE = 0.5),
      score_combine = "all",
      cache_dir = cache_dir
    ),
    tibble(
      marker_1 = "A", marker_2 = "C",
      uniprot_a = "P12345", uniprot_b = "P99999",
      ipSAE = 0.4, pDockQ2 = 0.5
    )
  )
})

test_that("extract_panel_interactions errors when scores are absent", {
  cache_dir <- tempfile("biogrid_noscore_")
  dir.create(cache_dir)
  raw_file <- file.path(
    tempfile("biogrid_"),
    "BIOGRID-ORGANISM-Homo_sapiens-5.0.259.tab3.txt"
  )
  dir.create(dirname(raw_file))
  utils::write.table(
    data.frame(
      `SWISS-PROT Accessions Interactor A` = "P12345",
      `SWISS-PROT Accessions Interactor B` = "Q67890",
      `Organism ID Interactor A` = 9606,
      `Organism ID Interactor B` = 9606,
      `Experimental System` = "Two-hybrid",
      `Experimental System Type` = "physical",
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    raw_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  build_biogrid_database(
    raw_file = raw_file,
    version = "5.0.259",
    cache_dir = cache_dir
  )
  expect_error(
    extract_panel_interactions(
      markers = c("A", "B"),
      database = "biogrid",
      marker_uniprot_map = marker_map,
      score_min = c(combined_score = 1),
      cache_dir = cache_dir
    ),
    "no score columns"
  )
})

test_that("create_marker_uniprot_map from Seurat assay", {
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

  expect_equal(
    create_marker_uniprot_map(seur, assay = "PNA") %>%
      arrange(marker, uniprot_id),
    tibble(
      marker = c("CD3e", "CD4", "CD4", "CD8"),
      uniprot_id = c("P07766", "P01730", "P01732", "P01732")
    )
  )
  expect_error(create_marker_uniprot_map(seur, assay = "missing"))
})
