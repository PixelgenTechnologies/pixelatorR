# Component ids are copied from a PXL file into Seurat cell names and then
# interpolated into DuckDB SQL. These tests require that text to be data:
# loading the file and running the public helpers must not run SQL that the
# component id contains. They fail while the names are glued or sprintf'd
# into the statement, and pass once each name is bound as a parameter.

benign_component <- "cellB"

write_crafted_pxl <- function(path, ids, edge_order = NULL, layout_order = NULL) {
  if (is.null(edge_order)) {
    edge_order <- ids
  }
  if (is.null(layout_order)) {
    layout_order <- ids
  }

  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE), path)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  DBI::dbExecute(con, "
    CREATE TABLE __adata__X (index VARCHAR, CD3 UINTEGER, CD4 UINTEGER)
  ")
  DBI::dbExecute(con, "
    CREATE TABLE __adata__obs (
      index VARCHAR,
      n_umi UBIGINT,
      n_edges UBIGINT,
      n_antibodies UBIGINT,
      isotype_fraction DOUBLE,
      average_k_core DOUBLE,
      reads_in_component BIGINT,
      n_umi1 UBIGINT,
      n_umi2 UBIGINT,
      sample VARCHAR
    )
  ")
  DBI::dbExecute(con, "
    CREATE TABLE __adata__var (
      index VARCHAR,
      antibody_count BIGINT,
      antibody_pct DOUBLE,
      components BIGINT,
      control VARCHAR,
      nuclear VARCHAR
    )
  ")
  DBI::dbExecute(con, "
    INSERT INTO __adata__var VALUES
      ('CD3', 3, 0.5, 2, 'no', 'no'),
      ('CD4', 3, 0.5, 2, 'no', 'no')
  ")
  DBI::dbExecute(con, "
    CREATE TABLE edgelist (
      marker_1 VARCHAR,
      marker_2 VARCHAR,
      umi1 UBIGINT,
      umi2 UBIGINT,
      read_count BIGINT,
      uei_count BIGINT,
      component VARCHAR
    )
  ")
  DBI::dbExecute(con, "
    CREATE TABLE layouts (
      index BIGINT,
      pixel_type VARCHAR,
      x DOUBLE,
      y DOUBLE,
      z DOUBLE,
      component VARCHAR,
      graph_projection VARCHAR,
      layout VARCHAR
    )
  ")
  DBI::dbExecute(con, "
    CREATE TABLE proximity (
      marker_1 VARCHAR,
      marker_2 VARCHAR,
      join_count UINTEGER,
      join_count_expected_mean DOUBLE,
      join_count_expected_sd DOUBLE,
      join_count_z DOUBLE,
      join_count_p DOUBLE,
      component VARCHAR
    )
  ")
  DBI::dbExecute(con, "CREATE TABLE metadata (value VARCHAR)")
  DBI::dbExecute(con, "
    INSERT INTO metadata VALUES
      ('{\"sample_name\":\"t\",\"version\":\"0.1.0\",\"technology\":\"single-cell-pna\"}')
  ")

  for (id in ids) {
    DBI::dbExecute(
      con,
      "INSERT INTO __adata__X VALUES (?, 2, 1)",
      params = list(id)
    )
    DBI::dbExecute(
      con,
      "INSERT INTO __adata__obs VALUES (?, 10, 2, 2, 0, 1, 10, 5, 5, 's')",
      params = list(id)
    )
    DBI::dbExecute(
      con,
      "INSERT INTO proximity VALUES ('CD3', 'CD4', 1, 1, 0, 0, 1, ?)",
      params = list(id)
    )
  }
  for (j in seq_along(edge_order)) {
    DBI::dbExecute(
      con,
      "INSERT INTO edgelist VALUES ('CD3', 'CD4', ?, ?, 5, 1, ?)",
      params = list(j * 10, j * 10 + 1, edge_order[[j]])
    )
  }
  for (j in seq_along(layout_order)) {
    DBI::dbExecute(
      con,
      "INSERT INTO layouts VALUES (?, 'A', 0.1, 0.2, 0.3, ?, 'bipartite', 'wpmds_3d')",
      params = list(j, layout_order[[j]])
    )
  }

  invisible(path)
}

# Run expr, then demand that it neither errored nor created `out`.
# A file created before a later error is still a failure.
expect_query_is_data_only <- function(out, expr) {
  if (file.exists(out)) {
    unlink(out)
  }
  err <- tryCatch(
    {
      force(expr)
      NULL
    },
    error = function(e) e
  )
  expect_false(file.exists(out))
  expect_null(err)
}

test_that("fetch_node_neighborhoods does not run SQL stored in a component id", {
  options(pixelatorR.verbose = FALSE)
  out <- tempfile(fileext = ".txt")
  on.exit(unlink(out), add = TRUE)
  payload <- paste0(
    "x') SELECT 1; COPY (SELECT 'NEIGH_SENTINEL') TO '", out,
    "'; WITH RECURSIVE comp_edges AS (",
    "SELECT umi1, umi2, marker_1, marker_2 FROM combined_edgelist ",
    "WHERE component = '", benign_component
  )
  pxl <- tempfile(fileext = ".pxl")
  on.exit(unlink(pxl), add = TRUE)
  write_crafted_pxl(pxl, c(payload, benign_component))

  se <- ReadPNA_Seurat(pxl, verbose = FALSE)
  expect_identical(colnames(se)[[1]], payload)

  expect_query_is_data_only(out, {
    fetch_node_neighborhoods(se, nodes_per_component = 1L, k = 1L)
  })
})

test_that("fetch_node_neighborhoods does not read local files through a component id", {
  options(pixelatorR.verbose = FALSE)
  secret <- tempfile(fileext = ".txt")
  writeLines("READ_SENTINEL_91", secret)
  on.exit(unlink(secret), add = TRUE)
  payload <- paste0(
    "x' OR CAST((SELECT content FROM read_text('", secret,
    "') LIMIT 1) AS INTEGER) = 1 OR component = '", benign_component
  )
  pxl <- tempfile(fileext = ".pxl")
  on.exit(unlink(pxl), add = TRUE)
  write_crafted_pxl(pxl, c(payload, benign_component))

  se <- ReadPNA_Seurat(pxl, verbose = FALSE)
  expect_identical(colnames(se)[[1]], payload)
  expect_no_error(
    fetch_node_neighborhoods(se, nodes_per_component = 1L, k = 1L)
  )
})

test_that("ReadPNA_layouts does not run SQL stored in a component id", {
  options(pixelatorR.verbose = FALSE)
  out <- tempfile(fileext = ".txt")
  on.exit(unlink(out), add = TRUE)
  # Stored after cellB, one row each, so the sorted-table branch is used.
  payload <- paste0(
    "z'); COPY (SELECT 'LAYOUT_SORTED') TO '", out,
    "'; SELECT CAST(index AS STRING) AS name, x, y, z FROM layouts ",
    "WHERE component IN ('", benign_component
  )
  pxl <- tempfile(fileext = ".pxl")
  on.exit(unlink(pxl), add = TRUE)
  write_crafted_pxl(pxl, c(benign_component, payload))

  expect_query_is_data_only(out, ReadPNA_layouts(pxl, verbose = FALSE))
})

test_that("ReadPNA_layouts does not run SQL from a component id when row ranges overlap", {
  options(pixelatorR.verbose = FALSE)
  out <- tempfile(fileext = ".txt")
  on.exit(unlink(out), add = TRUE)
  payload <- paste0(
    "z'; COPY (SELECT 'LAYOUT_NC') TO '", out,
    "'; SELECT CAST(index AS STRING) AS name, x, y, z FROM layouts ",
    "WHERE component = '", benign_component
  )
  pxl <- tempfile(fileext = ".pxl")
  on.exit(unlink(pxl), add = TRUE)
  write_crafted_pxl(
    pxl,
    c(benign_component, payload),
    layout_order = c(payload, benign_component, payload)
  )

  expect_query_is_data_only(out, ReadPNA_layouts(pxl, verbose = FALSE))
})

test_that("LoadCellGraphs does not run SQL stored in a component id", {
  options(pixelatorR.verbose = FALSE)
  out <- tempfile(fileext = ".txt")
  on.exit(unlink(out), add = TRUE)
  payload <- paste0(
    "z'; COPY (SELECT 'GRAPH_SORTED') TO '", out,
    "'; SELECT umi1, umi2, marker_1, marker_2 FROM edgelist ",
    "WHERE component = '", benign_component
  )
  pxl <- tempfile(fileext = ".pxl")
  on.exit(unlink(pxl), add = TRUE)
  write_crafted_pxl(pxl, c(benign_component, payload))
  se <- ReadPNA_Seurat(pxl, verbose = FALSE)

  expect_query_is_data_only(out, LoadCellGraphs(se, verbose = FALSE))
})

test_that("LoadCellGraphs does not run SQL from a component id when edges are interleaved", {
  options(pixelatorR.verbose = FALSE)
  out <- tempfile(fileext = ".txt")
  on.exit(unlink(out), add = TRUE)
  payload <- paste0(
    "z'; COPY (SELECT 'GRAPH_NC') TO '", out,
    "'; SELECT umi1, umi2, marker_1, marker_2 FROM edgelist ",
    "WHERE component = '", benign_component
  )
  pxl <- tempfile(fileext = ".pxl")
  on.exit(unlink(pxl), add = TRUE)
  write_crafted_pxl(
    pxl,
    c(benign_component, payload),
    edge_order = c(payload, benign_component, payload)
  )
  se <- ReadPNA_Seurat(pxl, verbose = FALSE)

  expect_query_is_data_only(out, LoadCellGraphs(se, verbose = FALSE))
})

test_that("components_marker_counts does not run SQL stored in a component id", {
  options(pixelatorR.verbose = FALSE)
  out <- tempfile(fileext = ".txt")
  on.exit(unlink(out), add = TRUE)
  payload <- paste0(
    "z')); COPY (SELECT 'MARKER_SORTED') TO '", out,
    "'; WITH src AS (SELECT umi1, umi2, marker_1, marker_2 FROM edgelist ",
    "WHERE component IN ('", benign_component
  )
  pxl <- tempfile(fileext = ".pxl")
  on.exit(unlink(pxl), add = TRUE)
  write_crafted_pxl(pxl, c(benign_component, payload))
  db <- PixelDB$new(pxl)
  on.exit(db$close(), add = TRUE)

  expect_query_is_data_only(out, db$components_marker_counts(NULL))
})

test_that("components_marker_counts does not run SQL from a component id when edges are interleaved", {
  options(pixelatorR.verbose = FALSE)
  out <- tempfile(fileext = ".txt")
  on.exit(unlink(out), add = TRUE)
  payload <- paste0(
    "z'); COPY (SELECT 'MARKER_NC') TO '", out,
    "'; WITH src AS (SELECT umi1, umi2, marker_1, marker_2 FROM edgelist ",
    "WHERE component = '", benign_component
  )
  pxl <- tempfile(fileext = ".pxl")
  on.exit(unlink(pxl), add = TRUE)
  write_crafted_pxl(
    pxl,
    c(benign_component, payload),
    edge_order = c(payload, benign_component, payload)
  )
  db <- PixelDB$new(pxl)
  on.exit(db$close(), add = TRUE)

  expect_query_is_data_only(out, db$components_marker_counts(NULL))
})
