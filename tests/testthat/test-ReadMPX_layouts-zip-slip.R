# Security spec for layout extraction.
#
# These tests fail while ReadMPX_layouts passes archive entry names through to
# utils::unzip. A fix must reject dangerous names before anything is written.
# LoadCellGraphs.MPXAssay is not called here. Its layout path is ReadMPX_layouts()
# with cells set to the component ids, which the selected-component test covers.
# ReadPNA_layouts reads layouts from the PXL database and does not unzip entry names.

write_u16 <- function(con, n) {
  writeBin(as.integer(n), con, size = 2L, endian = "little")
}

write_u32 <- function(con, n) {
  n <- as.numeric(n) %% 2^32
  writeBin(
    as.raw(c(
      n %% 256,
      floor(n / 256) %% 256,
      floor(n / 256^2) %% 256,
      floor(n / 256^3) %% 256
    )),
    con
  )
}

# Stored-method ZIP. R's internal unzip writes the local payload even when the
# CRC-32 field is zero, which is enough to observe where the bytes land.
write_stored_zip <- function(path, name, payload) {
  payload <- charToRaw(payload)
  name_n <- nchar(name, type = "bytes")
  con <- file(path, "wb")
  on.exit(close(con), add = TRUE)
  write_u32(con, 0x04034b50)
  write_u16(con, 20)
  write_u16(con, 0)
  write_u16(con, 0)
  write_u16(con, 0)
  write_u16(con, 33) # 1980-01-01
  write_u32(con, 0)
  write_u32(con, length(payload))
  write_u32(con, length(payload))
  write_u16(con, name_n)
  write_u16(con, 0)
  writeBin(charToRaw(name), con)
  writeBin(payload, con)
  cd_offset <- 30 + name_n + length(payload)
  write_u32(con, 0x02014b50)
  write_u16(con, 20)
  write_u16(con, 20)
  write_u16(con, 0)
  write_u16(con, 0)
  write_u16(con, 0)
  write_u16(con, 33)
  write_u32(con, 0)
  write_u32(con, length(payload))
  write_u32(con, length(payload))
  write_u16(con, name_n)
  write_u16(con, 0)
  write_u16(con, 0)
  write_u16(con, 0)
  write_u16(con, 0)
  write_u32(con, 0)
  write_u32(con, 0)
  writeBin(charToRaw(name), con)
  cd_size <- 46 + name_n
  write_u32(con, 0x06054b50)
  write_u16(con, 0)
  write_u16(con, 0)
  write_u16(con, 1)
  write_u16(con, 1)
  write_u32(con, cd_size)
  write_u32(con, cd_offset)
  write_u16(con, 0)
}

# R 4.5.1 and later drop ".." segments inside utils::unzip, but they still write
# the payload under the extraction directory. The temp-directory search is what
# fails on those versions until the package rejects the entry first.
find_payload <- function(root, payload) {
  if (!dir.exists(root)) {
    return(character())
  }
  files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE)
  files <- files[!dir.exists(files)]
  raw_payload <- charToRaw(payload)
  hits <- files[vapply(files, function(f) {
    info <- file.info(f)
    if (is.na(info$size) || info$size != length(raw_payload)) {
      return(FALSE)
    }
    identical(readBin(f, raw(), n = info$size), raw_payload)
  }, logical(1))]
  hits
}

# Parent-directory segments past the hive prefix and the session temp directory.
climb_out <- function(prefix, marker_rel) {
  paste0(prefix, paste(rep("..", 16), collapse = "/"), "/", marker_rel)
}

new_case <- function(label) {
  zip_dir <- "/tmp/pxl-zip-slip-zips"
  outside_root <- "/tmp/pxl-zip-slip-outside"
  dir.create(zip_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(outside_root, showWarnings = FALSE, recursive = TRUE)
  id <- paste0(label, "-", as.integer(Sys.time()), "-", sample.int(1e6, 1))
  marker <- file.path(outside_root, id, "marker")
  list(
    id = id,
    zip = file.path(zip_dir, paste0(id, ".pxl")),
    marker = marker,
    marker_rel = paste0("tmp/pxl-zip-slip-outside/", id, "/marker"),
    outside_root = outside_root,
    escape_dir = file.path(outside_root, id)
  )
}

expect_entry_not_materialized <- function(payload, case) {
  expect_false(file.exists(case$marker))
  expect_false(dir.exists(case$escape_dir))
  expect_length(find_payload(fs::path_temp(), payload), 0)
  expect_length(find_payload(tempdir(), payload), 0)
  expect_length(find_payload(case$outside_root, payload), 0)
}

prepare_archive <- function(case, entry, payload) {
  expect_lt(nchar(entry, type = "bytes"), 256)
  write_stored_zip(case$zip, entry, payload)
  listed <- utils::unzip(case$zip, list = TRUE)$Name
  expect_identical(listed, entry)
}

hive_prefix <- "layouts.parquet/graph_projection=bipartite/layout=pmds/component=CELL1/part-0/"

test_that("ReadMPX_layouts does not write outside when cells is NULL", {
  case <- new_case("null-cells")
  on.exit(unlink(c(case$zip, case$escape_dir, file.path(fs::path_temp(), "temp_layouts")), recursive = TRUE), add = TRUE)
  payload <- paste0("SENTINEL_NULL_", case$id)
  entry <- climb_out(hive_prefix, case$marker_rel)
  prepare_archive(case, entry, payload)

  expect_error(ReadMPX_layouts(case$zip, cells = NULL, verbose = FALSE))
  expect_entry_not_materialized(payload, case)
})

test_that("ReadMPX_layouts does not write outside for a selected component", {
  case <- new_case("selected-cell")
  on.exit(unlink(c(case$zip, case$escape_dir, file.path(fs::path_temp(), "temp_layouts")), recursive = TRUE), add = TRUE)
  payload <- paste0("SENTINEL_CELL_", case$id)
  entry <- climb_out(hive_prefix, case$marker_rel)
  prepare_archive(case, entry, payload)

  expect_error(ReadMPX_layouts(case$zip, cells = "CELL1", verbose = FALSE))
  expect_entry_not_materialized(payload, case)
})

test_that("ReadMPX_layouts does not write outside when the component name is missing", {
  case <- new_case("na-component")
  on.exit(unlink(c(case$zip, case$escape_dir, file.path(fs::path_temp(), "temp_layouts")), recursive = TRUE), add = TRUE)
  payload <- paste0("SENTINEL_NA_", case$id)
  entry <- climb_out("layouts.parquet/nohive/", case$marker_rel)
  prepare_archive(case, entry, payload)

  expect_error(ReadMPX_layouts(case$zip, cells = NULL, verbose = FALSE))
  expect_entry_not_materialized(payload, case)
})

test_that("ReadMPX_layouts does not extract an absolute layout entry name", {
  case <- new_case("absolute")
  on.exit(unlink(c(case$zip, case$escape_dir, file.path(fs::path_temp(), "temp_layouts")), recursive = TRUE), add = TRUE)
  payload <- paste0("SENTINEL_ABS_", case$id)
  entry <- paste0("/", case$marker_rel)
  # The absolute name still has to be chosen by the layouts.parquet substring search.
  entry <- sub("marker$", "layouts.parquet/marker", entry)
  case$marker <- sub("marker$", "layouts.parquet/marker", case$marker)
  prepare_archive(case, entry, payload)

  expect_error(ReadMPX_layouts(case$zip, cells = NULL, verbose = FALSE))
  expect_false(file.exists(case$marker))
  expect_length(find_payload(fs::path_temp(), payload), 0)
  expect_length(find_payload(tempdir(), payload), 0)
})

test_that("ReadMPX_layouts does not extract a layout entry name containing backslashes", {
  case <- new_case("backslash")
  on.exit(unlink(c(case$zip, case$escape_dir, file.path(fs::path_temp(), "temp_layouts")), recursive = TRUE), add = TRUE)
  payload <- paste0("SENTINEL_BACKSLASH_", case$id)
  entry <- paste0(
    "layouts.parquet\\",
    paste(rep("..", 8), collapse = "\\"),
    "\\",
    gsub("/", "\\", case$marker_rel, fixed = TRUE)
  )
  prepare_archive(case, entry, payload)

  expect_error(ReadMPX_layouts(case$zip, cells = NULL, verbose = FALSE))
  expect_entry_not_materialized(payload, case)
})

test_that("ReadMPX_layouts does not extract layouts for a component that was not requested", {
  case <- new_case("other-cell")
  on.exit(unlink(c(case$zip, case$escape_dir, file.path(fs::path_temp(), "temp_layouts")), recursive = TRUE), add = TRUE)
  payload <- paste0("SENTINEL_OTHER_", case$id)
  entry <- climb_out(hive_prefix, case$marker_rel)
  prepare_archive(case, entry, payload)

  expect_error(ReadMPX_layouts(case$zip, cells = "OTHER", verbose = FALSE))
  expect_entry_not_materialized(payload, case)
})

test_that("ReadMPX_layouts still extracts a relative layout entry", {
  case <- new_case("relative")
  on.exit(unlink(c(case$zip, file.path(fs::path_temp(), "temp_layouts")), recursive = TRUE), add = TRUE)
  payload <- paste0("SENTINEL_RELATIVE_", case$id)
  entry <- "layouts.parquet/graph_projection=bipartite/layout=pmds/component=CELL1/part-0.parquet"
  prepare_archive(case, entry, payload)

  tryCatch(ReadMPX_layouts(case$zip, cells = "CELL1", verbose = FALSE), error = function(e) e)
  hits <- unique(c(
    find_payload(fs::path_temp(), payload),
    find_payload(tempdir(), payload)
  ))
  expect_length(hits, 1)
  expect_false(file.exists(case$marker))
  expect_true(all(startsWith(normalizePath(hits), normalizePath(fs::path_temp()))))
})
