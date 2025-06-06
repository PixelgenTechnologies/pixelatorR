#' @include generics.R
NULL

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load methods
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param cells A character vector of cell names to load CellGraphs for
#' @param load_as Choose how the cell graph should be represented. This option has no
#' effect when \code{data_type = "PNA"} (see details below)
#' @param data_type One of "PNA" or "MPX"
#' @param add_marker_counts Should marker counts be added to the CellGraph objects?
#' Note that this option is ignored for the \code{data.frame} method when \code{data_type = "PNA"}.
#'
#' @param chunk_size Length of chunks used to load CellGraphs from edge list.
#' @param verbose Print messages
#'
#' @rdname LoadCellGraphs
#' @method LoadCellGraphs FileSystemDataset
#'
#' @export
#'
LoadCellGraphs.FileSystemDataset <- function(
  object,
  cells,
  load_as = c("bipartite", "Anode", "linegraph"),
  add_marker_counts = TRUE,
  data_type = c("PNA", "MPX"),
  chunk_size = 10,
  verbose = TRUE,
  ...
) {
  # Validate input parameters
  assert_vector(cells, type = "character", n = 1)
  assert_single_value(chunk_size, type = "integer")
  data_type <- match.arg(data_type, choices = c("PNA", "MPX"))

  # TODO: implement PNA method
  if (data_type == "PNA") {
    cli::cli_abort(
      c(
        "x" = "This method is not yet implemented for PNA data."
      )
    )
  }

  # Make sure that cells doesn't contain duplicated values
  if (sum(duplicated(cells)) > 0) {
    abort("'cells' cannot contain duplicated values.")
  }

  # Validate load_as
  load_as <- match.arg(load_as, choices = c("bipartite", "Anode", "linegraph"))

  # Select load function
  graph_load_fkn <-
    switch(load_as,
      "bipartite" = .load_mpx_as_bipartite,
      "Anode" = .load_mpx_as_anode,
      "linegraph" = .load_mpx_as_linegraph
    )

  # Convert edgelist to list of Cell Graphs
  if (verbose && check_global_verbosity()) {
    cli_alert("  Loading {length(cells)} edgelist(s) as {col_br_magenta(load_as)} graph(s)")
  }

  # Group cells ids into chunks
  cells_split <- split(cells, ceiling(seq_along(cells) / chunk_size))

  # Process chunks
  cellgraphs <- lapply(seq_along(cells_split), function(i) {
    cell_ids <- cells_split[[i]]

    # Load chunks for specific sample
    g_list <- try(
      {
        graph_load_fkn(object,
          cell_ids = cell_ids,
          add_markers = add_marker_counts
        )
      },
      silent = TRUE
    )

    if (inherits(g_list, what = "try-error") || any(sapply(g_list, is.null))) {
      cli::cli_abort(
        c("x" = "Failed to load edge list data. Most likely reason is that invalid cells were provided.")
      )
    }

    # Add marker counts
    if (add_marker_counts) {
      g_list <- lapply(g_list, function(g) {
        return(CreateCellGraphObject(g$graph, counts = g$counts, verbose = FALSE))
      })
    } else {
      g_list <- lapply(g_list, function(g) {
        return(CreateCellGraphObject(g, verbose = FALSE))
      })
    }

    return(g_list)
  }) %>% Reduce(c, .)

  return(cellgraphs)
}

#' Internal load method for MPX data
#'
#' @noRd
.load_mpx_graphs_from_df <- function(
  object,
  cells,
  load_as = c("bipartite", "Anode", "linegraph"),
  add_marker_counts = TRUE,
  verbose = TRUE
) {
  # Select load function
  graph_load_fkn <- switch(load_as,
    "bipartite" = .load_mpx_as_bipartite,
    "Anode" = .load_mpx_as_anode,
    "linegraph" = .load_mpx_as_linegraph
  )

  # Load cell graphs
  cg_list <- graph_load_fkn(object, cell_ids = cells, add_markers = add_marker_counts)

  if (add_marker_counts) {
    cg_list <- lapply(cg_list, function(g) {
      return(CreateCellGraphObject(g$graph, counts = g$counts, verbose = FALSE))
    })
  } else {
    cg_list <- lapply(cg_list, function(g) {
      return(CreateCellGraphObject(g, verbose = FALSE))
    })
  }

  return(cg_list)
}

#' Internal load method for PNA data
#'
#' @noRd
.load_pna_graphs_from_df <- function(
  object,
  cells,
  verbose = TRUE
) {
  # Select load function
  graph_load_fkn <- .load_pna_as_bipartite

  # Convert edgelist to list of Cell Graphs
  if (verbose && check_global_verbosity()) {
    cli_alert("  Loading {length(cells)} edgelist(s) as {col_br_magenta('bipartite')} graph(s)")
  }

  # Load chunks for specific sample
  g_list <- try(
    {
      graph_load_fkn(
        object,
        cell_ids = cells
      )
    },
    silent = TRUE
  )

  if (inherits(g_list, what = "try-error") || any(sapply(g_list, is.null))) {
    abort(glue("Failed to load edge list data. Most likely reason is that invalid cells were provided."))
  }

  # Create CellGraph objects
  g_list <- lapply(g_list, function(g) {
    return(CreateCellGraphObject(g, verbose = FALSE))
  })

  return(g_list)
}


#' @param add_layouts Load layouts from the PXL file if available.
#' @param force Force load graph(s) if they are already loaded
#' @param cl A cluster object created by makeCluster, or an integer
#' to indicate number of child-processes (integer values are ignored
#' on Windows) for parallel evaluations. See Details on performance
#' in the documentation for \code{pbapply}. The default is NULL,
#' which means that no parallelization is used.
#'
#' @rdname LoadCellGraphs
#' @method LoadCellGraphs MPXAssay
#'
#' @export
#'
LoadCellGraphs.MPXAssay <- function(
  object,
  cells = colnames(object),
  load_as = c("bipartite", "Anode", "linegraph"),
  add_marker_counts = TRUE,
  add_layouts = FALSE,
  force = FALSE,
  chunk_size = 10,
  cl = NULL,
  verbose = TRUE,
  ...
) {
  # Validate input parameters
  assert_vector(cells, type = "character", n = 1)

  # Make sure that cells doesn't contain duplicated values
  if (sum(duplicated(cells)) > 0) {
    cli::cli_abort(
      c(
        "i" = "{.var cells} cannot contain duplicated values.",
        "x" = "Cell IDs {.val {cells[duplicated(cells)]}} are duplicated"
      )
    )
  }
  assert_x_in_y(cells, colnames(object))
  load_as <- match.arg(load_as, choices = c("bipartite", "Anode", "linegraph"))

  # Check if cells are already loaded
  loaded_graphs <- !sapply(slot(object, name = "cellgraphs")[cells], is.null)
  if ((sum(!loaded_graphs) == 0) && !force) {
    if (verbose && check_global_verbosity()) {
      cli_alert_info("All cells are already loaded. Returning object unmodified.")
    }
    return(object)
  } else {
    if (force) {
      slot(object, "cellgraphs")[cells] <- rep(list(NULL), length(cells)) %>% set_names(nm = cells)
      loaded_graphs <- !sapply(slot(object, name = "cellgraphs")[cells], is.null)
    }
    cells_to_load <- setdiff(cells, cells[loaded_graphs])
    if (verbose && check_global_verbosity() && (length(cells_to_load) < length(cells))) {
      cli_alert_info(glue(
        "{length(cells) - length(cells_to_load)} CellGraphs loaded.",
        " Loading remaining {length(cells_to_load)} CellGraphs."
      ))
    }
  }

  # Find where components are located
  fs_map_unnested_filtered <- slot(object, name = "fs_map") %>%
    tidyr::unnest(cols = "id_map") %>%
    filter(current_id %in% cells)
  fs_map_nested_filtered <- fs_map_unnested_filtered %>%
    tidyr::nest(data = c("current_id", "original_id"))

  # Check for layouts
  if (add_layouts) {
    if (load_as != "bipartite") {
      cli::cli_abort(
        c("x" = "This function currently only supports bipartite graphs.")
      )
    }
    for (f in fs_map_nested_filtered$pxl_file) {
      pxl_file_info <- inspect_pxl_file(f)
      # Check if the file contains layouts
      if (!"layouts.parquet" %in% pxl_file_info$file_type) {
        cli::cli_abort(
          c("x" = "File {.file {f}} does not contain any component layouts.")
        )
      }
    }

    # Load precomputed layouts
    precomputed_layouts_sample <- list()
    layout_types_sample <- list()
    for (i in seq_len(nrow(fs_map_nested_filtered))) {
      f <- fs_map_nested_filtered$pxl_file[i]
      cell_id_data <- fs_map_nested_filtered$data[[i]] %>%
        filter(current_id %in% cells)
      original_id <- cell_id_data$original_id
      current_id <- cell_id_data$current_id
      precomputed_layouts <- ReadMPX_layouts(f, cells = original_id, graph_projection = load_as)
      precomputed_layouts <- lapply(precomputed_layouts, function(x) {
        x %>% set_names(nm = current_id)
      })
      layout_types_sample[[i]] <- names(precomputed_layouts)
      precomputed_layouts_sample[[i]] <- precomputed_layouts
    }

    # Keep union of all layout types
    all_layout_types <- Reduce(union, layout_types_sample)

    # Merge layout lists
    precomputed_layouts_merged <- list()
    for (layout_type in all_layout_types) {
      precomputed_layouts_merged[[layout_type]] <-
        Reduce(c, lapply(precomputed_layouts_sample, function(x) {
          x[[layout_type]]
        }))
    }
  }

  # Load cell graphs in chunks
  cg_list_full <- lapply(seq_len(nrow(fs_map_nested_filtered)), function(i) {
    f <- fs_map_nested_filtered$pxl_file[i]

    if (!fs::file_exists(f)) {
      cli::cli_abort(
        c(
          "x" = "File {.file {f}} does not exist.",
          "i" = "Run {.var ?RestorePaths} to get instructions on ",
          " " = "how to restore the PXL file paths."
        )
      )
    }

    # Unzip the edgelist parquet file to tmpdir
    utils::unzip(f, exdir = tempdir(), files = "edgelist.parquet")
    unz_pq_file <- file.path(tempdir(), "edgelist.parquet")
    pq_file <- fs::file_temp(ext = "parquet")

    # The file is renamed to make sure it has a unique path
    fs::file_move(unz_pq_file, pq_file)

    # Read the edgelist in memory, but only the necessary columns
    ar <- arrow::read_parquet(pq_file,
      col_select = c("upia", "upib", "marker", "component"),
      as_data_frame = FALSE
    )

    if (verbose && check_global_verbosity()) {
      cli_alert(glue(
        "   Loading CellGraphs for {nrow(fs_map_nested_filtered$data[[i]])} cells ",
        "from sample {fs_map_nested_filtered$sample[i]}"
      ))
    }

    # Split data into chunks determined by chunk_size
    id_data_chunks <- fs_map_nested_filtered$data[[i]] %>%
      mutate(group = ceiling(seq_len(n()) / chunk_size)) %>%
      rename(component = current_id) %>%
      group_by(group) %>%
      group_split()

    # Read cellgraphs
    cg_list <- pblapply(id_data_chunks, function(id_chunk) {
      # Filter edgelist data for the current chunk
      # Note that we use original_id which are the original
      # MPX component ids. current_id are the ids currenty used
      # for components in object
      edgelist_data <- ar %>%
        mutate(component = as.character(component)) %>%
        filter(component %in% id_chunk$original_id) %>%
        collect()

      # Send to internal data.frame method
      cg_list <- .load_mpx_graphs_from_df(
        edgelist_data,
        cells = id_chunk$original_id,
        load_as = load_as,
        add_marker_counts = add_marker_counts,
        verbose = verbose
      )

      return(cg_list)
    }, cl = cl) %>%
      unlist()

    # Update names of the cellgraph list
    cg_list <- cg_list[fs_map_nested_filtered$data[[i]]$original_id]
    cg_list <- set_names(cg_list, nm = fs_map_nested_filtered$data[[i]]$current_id)

    # Remove temporary file
    try_delete <- try(fs::file_delete(pq_file), silent = TRUE)
    if (inherits(try_delete, what = "try-error")) {
      cli_alert_warning("Failed to delete temporary edge list parquet file {pq_file}.")
    }
    if (inherits(try_delete, what = "try-error")) {
      cli_alert_warning("Failed to delete temporary edge list parquet file {pq_file}.")
    }

    return(cg_list)
  }) %>%
    unlist()

  # Add layouts to the list of cellgraphs if layouts were loaded
  if (add_layouts) {
    cg_list_full <- lapply(names(cg_list_full), function(nm) {
      cg <- cg_list_full[[nm]]
      cg@layout <- list()
      for (layout_type in all_layout_types) {
        if (attr(cg@cellgraph, "type") == "bipartite") {
          node_names <- rownames(cg@counts) %>%
            stringr::str_replace("-[A|B]", "")
        } else {
          node_names <- rownames(cg@counts)
        }
        # Rearrange layout node coordinates to match CellGraph node order
        coords <- precomputed_layouts_merged[[layout_type]][[nm]]
        coords <- coords[match(node_names, coords$name), ]
        cg@layout[[layout_type]] <- coords %>% select(-all_of("name"))
      }
      return(cg)
    }) %>%
      set_names(nm = names(cg_list_full))
  }

  # Fill cellgraphs slot list with the loaded CellGraphs
  slot(object, name = "cellgraphs")[fs_map_unnested_filtered$current_id] <- cg_list_full

  # Return object
  if (verbose && check_global_verbosity()) {
    cli_alert_success("Successfully loaded {length(cells_to_load)} CellGraph object(s).")
  }
  return(object)
}


#' @rdname LoadCellGraphs
#' @method LoadCellGraphs CellGraphAssay
#' @docType methods
#' @export
#'
LoadCellGraphs.CellGraphAssay <- LoadCellGraphs.MPXAssay

#' @rdname LoadCellGraphs
#' @method LoadCellGraphs CellGraphAssay5
#' @docType methods
#' @export
#'
LoadCellGraphs.CellGraphAssay5 <- LoadCellGraphs.MPXAssay


#' @rdname LoadCellGraphs
#' @method LoadCellGraphs PNAAssay
#'
#' @examples
#' # Read from PNAAssay
#' pxl_file <- minimal_pna_pxl_file()
#' seur_obj <- ReadPNA_Seurat(pxl_file)
#' pna_assay <- LoadCellGraphs(seur_obj[["PNA"]], cells = "0a45497c6bfbfb22")
#' CellGraphs(pna_assay)[["0a45497c6bfbfb22"]]
#'
#' @export
#'
LoadCellGraphs.PNAAssay <- function(
  object,
  cells = colnames(object),
  add_marker_counts = TRUE,
  add_layouts = FALSE,
  force = FALSE,
  chunk_size = 10,
  cl = NULL,
  verbose = TRUE,
  ...
) {
  # Validate input parameters
  assert_vector(cells, type = "character", n = 1)
  assert_unique(cells)
  assert_x_in_y(cells, colnames(object))

  # Check if cells are already loaded
  loaded_graphs <- !sapply(slot(object, name = "cellgraphs")[cells], is.null)
  if ((sum(!loaded_graphs) == 0) && !force) {
    if (verbose && check_global_verbosity()) {
      cli_alert_info("All cells are already loaded. Returning object unmodified.")
    }
    return(object)
  } else {
    if (force) {
      slot(object, "cellgraphs")[cells] <- rep(list(NULL), length(cells)) %>% set_names(nm = cells)
      loaded_graphs <- !sapply(slot(object, name = "cellgraphs")[cells], is.null)
    }
    cells_to_load <- setdiff(cells, cells[loaded_graphs])
    if (verbose && check_global_verbosity() && (length(cells_to_load) < length(cells))) {
      cli_alert_info(glue(
        "{length(cells) - length(cells_to_load)} CellGraphs loaded.",
        " Loading remaining {length(cells_to_load)} CellGraphs."
      ))
    }
  }

  # Find where components are located
  fs_map_unnested_filtered <- slot(object, name = "fs_map") %>%
    tidyr::unnest(cols = "id_map") %>%
    filter(current_id %in% cells)
  fs_map_nested_filtered <- fs_map_unnested_filtered %>%
    tidyr::nest(data = c("current_id", "original_id"))

  # Load cell graphs in chunks
  cg_list_full <- lapply(seq_len(nrow(fs_map_nested_filtered)), function(i) {
    f <- fs_map_nested_filtered$pxl_file[i]

    if (!fs::file_exists(f)) {
      abort(glue(
        "File '{col_br_blue(f)}' does not exist.\n",
        "Run ?RestorePaths to get instructions on ",
        "how to restore the PXL file paths."
      ))
    }

    # Split data into chunks determined by chunk_size
    id_data <- fs_map_nested_filtered$data[[i]] %>%
      mutate(group = ceiling(seq_len(n()) / chunk_size)) %>%
      rename(component = current_id)
    id_data_chunks <- split(id_data, id_data$group)

    # Setup connection to PXL database
    # TODO: can we speed this up with e.g. arrow?
    db <- PixelDB$new(f)
    on.exit(db$close())

    if (verbose && check_global_verbosity()) {
      cli_alert_info(
        "Fetching edgelists for {nrow(fs_map_nested_filtered$data[[i]])} cells ",
        "from sample {fs_map_nested_filtered$sample[i]}"
      )
    } else {
      opb <- getOption("pboptions")
      pbapply::pboptions(type = "none")
      on.exit(pbapply::pboptions(opb))
    }

    # Fetch the component edge tables from the data base
    edge_table_list <- pblapply(id_data_chunks, function(id_chunk) {
      db$components_edgelist(id_chunk$original_id, umi_data_type = "suffixed_string", include_all_columns = FALSE)
    })

    if (verbose && check_global_verbosity()) {
      cli_alert("Creating {.cls CellGraph} objects")
    }

    # Construct CellGraph objects from the edge tables
    cg_list <- pblapply(seq_along(edge_table_list), function(i) {
      edge_table <- edge_table_list[[i]]
      id_chunk <- id_data_chunks[[i]]

      # Send to internal data.frame method
      cg_list <- .load_pna_graphs_from_df(
        edge_table,
        cells = id_chunk$original_id,
        verbose = FALSE
      )

      return(cg_list)
    }, cl = cl) %>%
      unlist()

    # Load markers counts if specified
    if (add_marker_counts) {
      if (verbose && check_global_verbosity()) {
        cli_alert("Fetching marker counts")
      }
      # Fetch marker counts from the data base
      # The marker counts matrices are created from the edgelist
      # in the data base using SQL queries
      marker_counts_list <- pblapply(id_data_chunks, function(id_chunk) {
        db$components_marker_counts(id_chunk$original_id, as_sparse = TRUE)
      }) %>% Reduce(c, .)

      if (verbose && check_global_verbosity()) {
        cli_alert("Adding marker counts to {.cls CellGraph} object(s)")
      }
      # Fill marker counts slot list with the loaded marker counts
      cg_list <- pblapply(names(cg_list), function(nm) {
        cg <- cg_list[[nm]]
        counts <- marker_counts_list[[nm]]
        counts <- counts[match(cg@cellgraph %N>% pull(name), rownames(counts)), ]
        cg@counts <- counts
        return(cg)
      }) %>%
        set_names(nm = names(cg_list))
    }

    # Load layout if specified
    if (add_layouts) {
      if (verbose && check_global_verbosity()) {
        cli_alert("Fetching layouts")
      }
      if (!"layouts" %in% db$names()) {
        cli::cli_abort(c("x" = "Layouts are missing from the PXL file."))
      }
      layout_list <- pblapply(id_data_chunks, function(id_chunk) {
        db$components_layout(id_chunk$original_id, verbose = FALSE)
      }) %>% Reduce(c, .)

      if (verbose && check_global_verbosity()) {
        cli_alert("Adding layouts to {.cls CellGraph} object(s)")
      }

      cg_list <- pblapply(names(cg_list), function(nm) {
        cg <- cg_list[[nm]]
        layout <- layout_list[[nm]]
        layout <- layout[match(cg@cellgraph %N>% pull(name), layout$name), ] %>% select(-name)
        cg@layout <- list(wpmds_3d = layout)
        return(cg)
      }) %>%
        set_names(nm = names(cg_list))
    }

    # Update names of the cellgraph list
    cg_list <- cg_list[fs_map_nested_filtered$data[[i]]$original_id]
    cg_list <- set_names(cg_list, nm = fs_map_nested_filtered$data[[i]]$current_id)

    return(cg_list)
  }) %>%
    unlist()

  # Fill cellgraphs slot list with the loaded CellGraphs
  slot(object, name = "cellgraphs")[fs_map_unnested_filtered$current_id] <- cg_list_full

  # Return object
  if (verbose && check_global_verbosity()) {
    cli_alert_success("Successfully loaded {length(cells_to_load)} {.cls CellGraph} object(s).")
  }
  return(object)
}

#' @rdname LoadCellGraphs
#' @method LoadCellGraphs PNAAssay5
#' @docType methods
#' @export
#'
LoadCellGraphs.PNAAssay5 <- LoadCellGraphs.PNAAssay


#' @param assay Assay name
#' @param cells A character vector of cell names to load CellGraphs for
#' @param verbose Print messages
#'
#' @rdname LoadCellGraphs
#' @method LoadCellGraphs Seurat
#'
#' @export
#'
LoadCellGraphs.Seurat <- function(
  object,
  assay = NULL,
  cells = colnames(object),
  load_as = c("bipartite", "Anode", "linegraph"),
  add_marker_counts = TRUE,
  add_layouts = FALSE,
  force = FALSE,
  chunk_size = 10,
  cl = NULL,
  verbose = TRUE,
  ...
) {
  # Use default assay if assay = NULL
  if (!is.null(assay)) {
    assert_single_value(assay, type = "string")
  } else {
    # Use default assay if assay = NULL
    assay <- DefaultAssay(object)
  }

  # Validate assay
  cg_assay <- object[[assay]]
  if (!inherits(cg_assay, c("MPXAssay", "PNAAssay", "PNAAssay5"))) {
    cli::cli_abort(
      c(
        "x" = "Invalid assay type {.cls {class(cg_assay)}}. Expected one of:",
        " " = "{.cls {c('CellGraphAssay', 'CellGraphAssay5', 'PNAAssay', 'PNAAssay5')}}"
      )
    )
  }

  # Load cell graphs
  cg_assay <-
    LoadCellGraphs(
      cg_assay,
      cells = cells,
      load_as = load_as,
      add_marker_counts = add_marker_counts,
      add_layouts = add_layouts,
      force = force,
      chunk_size = chunk_size,
      cl = cl,
      verbose = verbose,
      ...
    )

  object[[assay]] <- cg_assay

  return(object)
}

#' Load bipartite graph from edgelist
#'
#' @param edge_table A \code{data.frame}
#' @param cell_ids An integer vector of cell IDs
#'
#' @noRd
.load_pna_as_bipartite <- function(
  edge_table,
  cell_ids
) {
  components <- pull(edge_table, component)

  component_graphs <- lapply(cell_ids, function(i) {
    edge_table_cur <- edge_table[components == i, ] %>% select(-component)

    # Convert to a tbl_graph
    g <- edge_table_cur %>%
      select(umi1, umi2) %>%
      as.matrix() %>%
      igraph::graph_from_edgelist(directed = FALSE)
    g <- as_tbl_graph(g, directed = FALSE) %>%
      mutate(node_type = stringr::str_extract(name, "umi[1|2]"))

    # Return results
    attr(g, "type") <- "bipartite"
    attr(g, "component_id") <- i
    attr(g, "assay_type") <- "PNA"
    return(g)
  }) %>%
    set_names(nm = cell_ids)

  return(component_graphs[cell_ids])
}


#' Load bipartite graph from edgelist
#'
#' @param arrow_data An arrow dataset
#' @param cell_ids A character vector of cell IDs
#'
#' @noRd
.load_mpx_as_bipartite <- function(
  arrow_data,
  cell_ids,
  add_markers = TRUE
) {
  edge_table <- arrow_data %>%
    mutate(component = as.character(component)) %>%
    filter(component %in% cell_ids) %>%
    collect() %>%
    group_by(component)

  edge_table_split <- edge_table %>%
    group_split() %>% # Split into list, with one element per component
    set_names(nm = edge_table %>% dplyr::group_data() %>% pull(component)) # Name list with component IDs

  edge_table_split <- lapply(names(edge_table_split), function(nm) {
    edge_table <- edge_table_split[[nm]]

    # Add a suffix to upia / upib and convert to a tbl_graph
    g <- edge_table %>%
      select(upia, upib, marker) %>%
      mutate(across(where(is.factor), as.character)) %>%
      mutate(upia = paste0(upia, "-A"), upib = paste0(upib, "-B")) %>%
      as_tbl_graph(directed = FALSE)

    # Add node type as node attribute
    g <- g %>%
      # Replace node attributes with name / node_type
      mutate(
        !!!g %>%
          as_tibble() %>%
          mutate(id = name) %>%
          tidyr::separate(col = name, into = c("name", "node_type"))
      ) %>%
      select(-name) %>%
      rename(name = id) %>%
      select(name, node_type)

    if (add_markers) {
      # Aggregate counts per upia / upib combination
      markers_s <- g %E>%
        as_tibble() %>%
        group_by(from, to, marker) %>%
        summarize(n = n(), .groups = "drop")

      # Count markers per upia/upib pair and fill empty cells with 0
      markers_wide <- markers_s %>%
        tidyr::pivot_wider(id_cols = c("from", "to"), names_from = marker, values_from = n, values_fill = 0)

      # Get new edges from pivot table (no duplicated edges left)
      new_edges <- markers_wide %>% select(from, to)

      # Get edge count matrix from pivot table
      edge_cntMatrix <- as(markers_wide[, 3:ncol(markers_wide)] %>% as.matrix(), "dgCMatrix")

      # Create node count matrix by aggregating counts for A nodes and B nodes
      grp1 <- markers_wide %>% pull(from)
      upia_markers <- as(markers_wide %>% select(-from, -to) %>% rowsum(group = grp1) %>% as.matrix(), "dgCMatrix")
      grp2 <- markers_wide %>% pull(to)
      upib_markers <- as(markers_wide %>% select(-from, -to) %>% rowsum(group = grp2) %>% as.matrix(), "dgCMatrix")
      node_cntMatrix <- rbind(upia_markers, upib_markers)
      rownames(node_cntMatrix) <- g %>% pull(name)

      # Replace with collapsed edges and return a simple undirected graph
      g <- g %E>% # Activate edges
        select(-marker) %>%
        filter(FALSE) %>% # Remove all edges
        bind_edges(new_edges) %N>% # Add the distinct edges matching the count matrix and reactivate nodes
        select(name, node_type)
    } else {
      # Skip marker counts and return a simple undirected graph
      new_edges <- g %E>% as_tibble() %>%
        select(from, to) %>%
        distinct()

      # Replace with collapsed edges and return a simple undirected graph
      g <- g %E>%
        select(-marker) %>%
        filter(FALSE) %>% # Remove all edges
        bind_edges(new_edges) %N>% # Add distinct edges and reactivate nodes
        select(name, node_type)
    }

    # Return results
    attr(g, "type") <- "bipartite"
    attr(g, "component_id") <- nm
    attr(g, "assay_type") <- "MPX"
    if (add_markers) {
      return(list(graph = g, counts = node_cntMatrix, counts_edges = edge_cntMatrix))
    } else {
      return(g)
    }
  }) %>%
    set_names(nm = names(edge_table_split))

  return(edge_table_split[cell_ids])
}


#' Load linegraph from edgelist
#'
#' @param arrow_data An arrow dataset
#' @param cell_id A character vector of cell IDs
#'
#' @noRd
.load_mpx_as_linegraph <- function(
  arrow_data,
  cell_ids,
  add_markers = TRUE
) {
  # Start by loading graph as bipartite
  g_list <- .load_mpx_as_bipartite(arrow_data, cell_ids = cell_ids, add_markers = add_markers)

  g_list <- lapply(names(g_list), function(nm) {
    g <- g_list[[nm]]

    # Unpack list elements if needed
    if (is.list(g)) {
      g_bipartite <- g$graph
    } else {
      g_bipartite <- g
    }

    # Convert squished bipartite graph to a linegraph
    g_line <- g_bipartite %>%
      tidygraph::to_linegraph()
    g_line <- g_line[[1]]

    # Return
    attr(g_line, "type") <- "linegraph"
    attr(g_line, "component_id") <- nm
    if (add_markers) {
      # For the linegraph, the node counts are the same as
      # the edge counts of the bipartite graph
      return(list(graph = g_line, counts = g$counts_edges))
    } else {
      return(g_line)
    }
  }) %>%
    set_names(nm = names(g_list))

  return(g_list)
}


#' Load A-node projected graph from edgelist
#'
#' @param arrow_data An arrow dataset
#' @param cell_id A cell ID
#'
#' @noRd
.load_mpx_as_anode <- function(
  arrow_data,
  cell_ids,
  add_markers = TRUE
) {
  # Fetch edgelist from parquet file and group by component
  edge_table <- arrow_data %>%
    mutate(component = as.character(component)) %>%
    filter(component %in% cell_ids) %>%
    collect() %>%
    group_by(component)

  # Split into list with one component per element
  edge_table_split <- edge_table %>%
    group_split() %>%
    set_names(nm = edge_table %>% dplyr::group_data() %>% pull(component))

  edge_table_split <- lapply(names(edge_table_split), function(nm) {
    edge_table <- edge_table_split[[nm]] %>%
      select(upia, upib, marker) %>%
      mutate(across(where(is.factor), as.character))

    # Anode projection
    g_anode <- edge_table %>%
      left_join(edge_table,
        by = c("upib"),
        relationship = "many-to-many",
        suffix = c("1", "2")
      ) %>%
      select(upia1, upia2) %>%
      distinct() %>%
      filter(upia1 < upia2)

    # Convert edge list to tbl_graph (undirected)
    g_anode <- g_anode %>%
      as_tbl_graph(directed = FALSE)

    if (add_markers) {
      # Group by upia and marker
      markers_s <- edge_table %>%
        as_tibble() %>%
        group_by(upia, marker) %>%
        summarize(n = n(), .groups = "drop")

      # Pivot marker counts into a wide format
      markers_wide <- markers_s %>%
        pivot_wider(id_cols = upia, names_from = marker, values_from = n, values_fill = 0)

      cntMatrix <- as(markers_wide[, 2:ncol(markers_wide)] %>% as.matrix(), "dgCMatrix")

      # Sort matrix rows to match graph names
      nds <- g_anode %>%
        pull(name)
      cntMatrix <- cntMatrix[match(nds, markers_wide$upia), ]
      rownames(cntMatrix) <- nds
    }

    # Return results
    attr(g_anode, "type") <- "Anode"
    attr(g_anode, "component_id") <- nm
    if (add_markers) {
      return(list(graph = g_anode, counts = cntMatrix))
    } else {
      return(g)
    }
  }) %>%
    set_names(nm = names(edge_table_split))

  # Return list with the correct order of components
  return(edge_table_split[cell_ids])
}
