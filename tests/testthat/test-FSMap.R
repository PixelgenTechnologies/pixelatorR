seur_obj_mpx <- ReadMPX_Seurat(minimal_mpx_pxl_file())
seur_obj_pna <- ReadPNA_Seurat(minimal_pna_pxl_file())

test_that("FSMap works as expected", {
  # Seurat object
  expect_no_error(fs_map <- FSMap(seur_obj_mpx))
  expect_equal(dim(fs_map), c(1, 3))
  expect_equal(names(fs_map), c("id_map", "sample", "pxl_file"))
  expect_no_error(FSMap(seur_obj_mpx) <- FSMap(seur_obj_mpx))

  # CellGraphAssay
  expect_no_error(fs_map_mpx <- FSMap(seur_obj_mpx[["mpxCells"]]))
  expect_equal(dim(fs_map_mpx), c(1, 3))
  expect_equal(names(fs_map_mpx), c("id_map", "sample", "pxl_file"))
  expect_no_error(FSMap(seur_obj_mpx[["mpxCells"]]) <- FSMap(seur_obj_mpx[["mpxCells"]]))

  # PNAAssay
  expect_no_error(fs_map_pna <- FSMap(seur_obj_pna[["PNA"]]))
  expect_equal(dim(fs_map_pna), c(1, 3))
  expect_equal(names(fs_map_pna), c("id_map", "sample", "pxl_file"))
  expect_no_error(FSMap(seur_obj_pna[["PNA"]]) <- FSMap(seur_obj_pna[["PNA"]]))

  # Check that FSMap can be updated
  tmp_path <- fs::file_temp(ext = "pxl")
  fs::file_copy(fs_map$pxl_file, tmp_path)
  expect_no_error({
    FSMap(seur_obj_mpx) <- FSMap(seur_obj_mpx) |>
      mutate(pxl_file = tmp_path |> as.character())
  })

  # Check that the update fails if the file is missing
  fs::file_delete(tmp_path)
  expect_error({
    FSMap(seur_obj_mpx) <- FSMap(seur_obj_mpx) |>
      mutate(pxl_file = tmp_path |> as.character())
  })

  # Check that FSMap works with multiple PNA assays
  expect_warning(seur_obj_pna[["PNA2"]] <- seur_obj_pna[["PNA"]])
  expect_warning(FSMap(seur_obj_pna), "Multiple PNA/MPX assays found")
  expect_equal(dim(FSMap(seur_obj_pna[["PNA"]])), c(1, 3))
  expect_equal(FSMap(seur_obj_pna[["PNA"]]), FSMap(seur_obj_pna[["PNA2"]]))
})

test_that("FSMap fails when invalid input is provided", {
  expect_error(fs_map <- FSMap("Invalid"))
  expect_error(FSMap(seur_obj) <- "Invalid")
})

test_that("FSMap S3 subclass and backend generics work and preserve class across PNAAssay operations", {
  for (assay_version in c("v3", "v5")) {
    options(Seurat.object.assay.version = assay_version)
    on.exit(options(Seurat.object.assay.version = "v3"), add = TRUE)

    seur1 <- suppressWarnings(ReadPNA_Seurat(minimal_pna_pxl_file(), verbose = FALSE))
    seur2 <- suppressWarnings(ReadPNA_Seurat(minimal_pna_pxl_file(), verbose = FALSE))

    custom_fs_map <- function(fs_map, uri = "mock://sample1.pxl") {
      fs_map$pxl_file <- uri
      class(fs_map) <- c("mock_fs_map", class(fs_map))
      fs_map
    }

    local_bindings(
      validate_fs_map_paths.mock_fs_map = function(fs_map, call = rlang::caller_env(), ...) {
        bad <- !startsWith(fs_map$pxl_file, "mock://")
        if (any(bad)) {
          cli::cli_abort("Invalid mock URI", call = call)
        }
        invisible(NULL)
      },
      lazy_load_table.mock_fs_map = function(fs_map,
                                             table_name = "proximity",
                                             calc_log2ratio = TRUE,
                                             proteins_keep = NULL,
                                             union = TRUE,
                                             call = rlang::caller_env(),
                                             ...) {
        real_map <- fs_map
        real_map$pxl_file <- rep(minimal_pna_pxl_file(), nrow(real_map))
        class(real_map) <- setdiff(class(real_map), "mock_fs_map")
        # Deduplicate temp copies if multiple rows
        if (nrow(real_map) > 1L) {
          tmp_files <- vapply(seq_len(nrow(real_map)), function(i) {
            tf <- fs::file_temp(ext = "pxl")
            fs::file_copy(minimal_pna_pxl_file(), tf, overwrite = TRUE)
            as.character(tf)
          }, character(1))
          real_map$pxl_file <- tmp_files
        }
        lazy_load_table.default(
          fs_map = real_map,
          table_name = table_name,
          calc_log2ratio = calc_log2ratio,
          proteins_keep = proteins_keep,
          union = union,
          call = call,
          ...
        )
      },
      fetch_cell_graphs.mock_fs_map = function(fs_map,
                                               cells_to_load,
                                               add_layouts = FALSE,
                                               verbose = TRUE,
                                               ...) {
        real_map <- fs_map
        real_map$pxl_file <- rep(minimal_pna_pxl_file(), nrow(real_map))
        class(real_map) <- setdiff(class(real_map), "mock_fs_map")
        fetch_cell_graphs.default(
          fs_map = real_map,
          cells_to_load = cells_to_load,
          add_layouts = add_layouts,
          verbose = verbose,
          ...
        )
      },
      .env = globalenv()
    )

    FSMap(seur1) <- custom_fs_map(FSMap(seur1), "mock://s1.pxl")
    FSMap(seur2) <- custom_fs_map(FSMap(seur2), "mock://s2.pxl")
    expect_s3_class(FSMap(seur1), "mock_fs_map")

    # Custom validator rejects non-mock URI
    bad_map <- FSMap(seur1)
    bad_map$pxl_file <- "/nonexistent/local.pxl"
    expect_error(FSMap(seur1) <- bad_map, "Invalid mock URI")

    # RenameCells preserves S3 class
    renamed_assay <- RenameCells(
      seur1[["PNA"]],
      new.names = paste0("S1_", colnames(seur1))
    )
    expect_s3_class(FSMap(renamed_assay), "mock_fs_map")

    # merge preserves S3 class when all inputs share it
    seur_merged <- merge(seur1, seur2, add.cell.ids = c("S1", "S2"))
    expect_s3_class(FSMap(seur_merged), "mock_fs_map")
    expect_equal(FSMap(seur_merged)$pxl_file, c("mock://s1.pxl", "mock://s2.pxl"))

    # subset preserves S3 class (both with cells and with features only)
    seur_sub_cells <- suppressWarnings(subset(seur_merged, cells = colnames(seur_merged)[1:3]))
    expect_s3_class(FSMap(seur_sub_cells), "mock_fs_map")
    expect_equal(nrow(FSMap(seur_sub_cells)), 1L)

    seur_sub_feats <- suppressWarnings(subset(seur_merged[["PNA"]], features = rownames(seur_merged)[1:10]))
    expect_s3_class(FSMap(seur_sub_feats), "mock_fs_map")
    expect_equal(nrow(FSMap(seur_sub_feats)), 2L)

    # Merging assays with incompatible fs_map backend classes aborts regardless of input order
    seur_local <- ReadPNA_Seurat(minimal_pna_pxl_file(), load_proximity_scores = FALSE, verbose = FALSE)
    expect_error(
      merge(seur1, seur_local, add.cell.ids = c("A", "B")),
      "incompatible `fs_map` backend classes"
    )
    expect_error(
      merge(seur_local, seur1, add.cell.ids = c("A", "B")),
      "incompatible `fs_map` backend classes"
    )

    # JoinLayers preserves S3 class on PNAAssay5
    if (assay_version == "v5") {
      joined_assay <- JoinLayers(seur_merged[["PNA"]])
      expect_s3_class(FSMap(joined_assay), "mock_fs_map")
    }

    # ProximityScores(lazy = TRUE) dispatches through lazy_load_table and join_proximity_metadata
    prox_lazy <- ProximityScores(
      seur_merged,
      lazy = TRUE,
      add_marker_proportions = TRUE,
      meta_data_columns = "n_umi"
    )
    expect_true("pxl_all_counts" %in% DBI::dbListTables(prox_lazy$src$con))
    prox_df <- collect(prox_lazy)
    DBI::dbDisconnect(prox_lazy$src$con)
    expect_true(all(c("count_1", "count_2", "p1", "p2", "n_umi") %in% colnames(prox_df)))
    expect_gt(nrow(prox_df), 0L)

    # Edgelists dispatches through lazy_load_table and join_proximity_metadata
    el_df <- Edgelists(seur_merged, lazy = FALSE, meta_data_columns = "n_umi")
    expect_true("n_umi" %in% colnames(el_df))
    expect_gt(nrow(el_df), 0L)

    # LoadCellGraphs dispatches through fetch_cell_graphs
    target_cell <- colnames(seur_merged)[1]
    seur_loaded <- LoadCellGraphs(seur_merged, cells = target_cell, verbose = FALSE)
    expect_s4_class(CellGraphs(seur_loaded)[[target_cell]], "CellGraph")
  }
})
