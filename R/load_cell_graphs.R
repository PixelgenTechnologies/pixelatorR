#' @include generics.R
NULL

# Declarations used in package check
globalVariables(
  names = c('g', 'from', 'to', 'node_type'),
  package = 'pixelatorR',
  add = TRUE
)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load methods
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @param cells A character vector of cell names to load CellGraphs for
#' @param load_as Choose how the cell graph should be represented (see details below)
#' @param add_marker_counts Should marker counts be added to the CellGraph objects?
#' @param force Force load graph(s) if they are already loaded
#' @param chunk_length Length of chunks used to load CellGraphs from edge list.
#' Smaller chunks will likely take longer time to load, but decreases memory usage.
#' @param verbose Print messages
#' @import cli
#' @import glue
#'
#' @rdname LoadCellGraphs
#' @method LoadCellGraphs CellGraphAssay
#'
#' @export
#'
LoadCellGraphs.CellGraphAssay <- function (
  object,
  cells = colnames(object),
  load_as = c("bipartite", "Anode", "linegraph"),
  add_marker_counts = TRUE,
  force = FALSE,
  chunk_length = 10,
  verbose = TRUE,
  ...
) {

  # Validate input parameters
  stopifnot(
    "'cells' must be a non-empty character vector of cell names" =
      is.character(cells) &&
      (length(cells) > 0)
  )
  stopifnot(
    "'cells' must be present in 'object'" =
      all(cells %in% colnames(object))
  )
  load_as <- match.arg(load_as, choices = c("bipartite", "Anode", "linegraph"))

  # Check if cells are already loaded
  loaded_graphs <- !sapply(slot(object, name = "cellgraphs")[cells], is.null)
  if ((sum(!loaded_graphs) == 0) & !force) {
    if (verbose && check_global_verbosity())
      cli_alert_info("All cells are alredy loaded. Returning CellGraphAssay unmodified.")
    return(object)
  } else {
    if (force) {
      slot(object, "cellgraphs")[cells] <- rep(list(NULL), length(cells)) %>% setNames(nm = cells)
      loaded_graphs <- !sapply(slot(object, name = "cellgraphs")[cells], is.null)
    }
    cells_to_load <- setdiff(cells, cells[loaded_graphs])
    if (verbose && check_global_verbosity() & (length(cells_to_load) < length(cells)))
      cli_alert_info(glue("{length(cells) - length(cells_to_load)} CellGraphs loaded.",
                          " Loading remaining {length(cells_to_load)} CellGraphs."))
  }

  # Check for valid arrow_dir
  arrow_dir <- slot(object, name = "arrow_dir")
  if (length(arrow_dir) > 1) {
    abort("'arrow_dir' should be a single directory")
  }
  if (is.na(arrow_dir)) {
    cli_alert_warning("'arrow_dir' is missing from object")
  } else {
    stopifnot("'arrow_dir' must be a directory" = dir.exists(arrow_dir))
  }

  # Check for valid arrow_data
  arrow_data <- slot(object, name = "arrow_data")
  msg <- tryCatch(arrow_data %>% nrow(), error = function(e) "Error")
  msg <- msg %||% "Error"
  if (msg == "Error") {
    cli_alert_warning("'arrow_data' is either missing from the object or dead. Attempting to load edgelist from stored 'arrow_dir'.")
    if (length(arrow_dir) == 0) {
      cli_alert_info("'arrow_dir' is missing from object. Returning object unmodified.")
    }
    arrow_data <- ReadMPX_arrow_edgelist(path = arrow_dir, outdir = arrow_dir, verbose = FALSE)
  }

  # Ensure that all cell names are available in edgelist
  stopifnot(
    "All 'cells' must be present in the edgelist" =
      all(cells_to_load %in% (arrow_data %>% pull(component, as_vector = TRUE)))
  )

  graph_load_fkn <- switch (load_as,
    "bipartite" = .load_as_bipartite,
    "Anode" = .load_as_anode,
    "linegraph" = .load_as_linegraph
  )

  # Convert edgelist to list of Cell Graphs
  if (verbose && check_global_verbosity())
    cli_alert("  Loading {length(cells_to_load)} edgelist(s) as {col_br_magenta(load_as)} graph(s)")

  # Split cells id into chunks
  cells_to_load_chunks <- split(cells_to_load, ceiling(seq_along(cells_to_load) / chunk_length))

  # Load cell graphs
  p <- progressr::progressor(along = cells_to_load_chunks)
  cellgraphs <- lapply(cells_to_load_chunks, function (cell_ids) {
    g_list <- graph_load_fkn(arrow_data, cell_ids = cell_ids, add_markers = add_marker_counts)
    if (add_marker_counts) {
      g_list <- lapply(g_list, function(g) {
        return(CreateCellGraphObject(g$graph, counts = g$counts, verbose = FALSE))
      })
    } else {
      g_list <- lapply(g_list, function(g) {
        return(CreateCellGraphObject(g, verbose = FALSE))
      })
    }
    p()
    return(g_list)
  }) %>% Reduce(c, .)

  # Fill empty elements if needed
  slot(object, name = "cellgraphs")[cells_to_load] <- cellgraphs[cells_to_load]

  # Return object
  if (verbose && check_global_verbosity())
    cli_alert_success("Successfully loaded {length(cells_to_load)} CellGraph object(s).")
  return(object)
}


#' @param assay Assay name
#' @param cells A character vector of cell names to load CellGraphs for
#' @param verbose Print messages
#'
#' @rdname LoadCellGraphs
#' @method LoadCellGraphs Seurat
#'
#' @export
#'
LoadCellGraphs.Seurat <- function (
  object,
  assay = NULL,
  cells = colnames(object),
  load_as = c("bipartite", "Anode", "linegraph"),
  add_marker_counts = TRUE,
  force = FALSE,
  chunk_length = 10,
  verbose = TRUE,
  ...
) {

  # Use default assay if assay = NULL
  if (!is.null(assay)) {
    stopifnot("'assay' must be a character of length 1" = is.character(assay) & (length(assay) == 1))
  } else {
    # Use default assay if assay = NULL
    assay <- DefaultAssay(object)
  }

  # Validate assay
  cg_assay <- object[[assay]]
  if (!inherits(cg_assay, what = "CellGraphAssay")) {
    abort(glue("Invalid assay type '{class(cg_assay)}'. Expected a 'CellGraphAssay'"))
  }

  # Load cell graphs
  cg_assay <-
    LoadCellGraphs(
      cg_assay,
      cells = cells,
      load_as = load_as,
      add_marker_counts = add_marker_counts,
      force = force,
      chunk_length = chunk_length,
      verbose = verbose,
      ...
    )

  object[[assay]] <- cg_assay

  return(object)
}


#' Load bipartite graph from edgelist
#'
#' @param arrow_data An arrow dataset
#' @param cell_ids A character vector of cell IDs
#'
#' @import dplyr
#' @importFrom tidygraph as_tbl_graph `%E>%` `%N>%`
#'
#' @noRd
.load_as_bipartite <- function (
  arrow_data,
  cell_ids,
  add_markers = TRUE
) {
  edge_table <- arrow_data %>%
    filter(component %in% cell_ids) %>%
    collect() %>%
    group_by(component)

  edge_table_split <- edge_table %>%
    group_split() %>% # Split into list, with one element per component
    setNames(nm = edge_table %>% dplyr::group_data() %>% pull(component)) # Name list with component IDs

  edge_table_split <- lapply(edge_table_split, function(edge_table) {

    # Add a suffix to upia / upib and convert to a tbl_graph
    g <- edge_table %>%
      select(upia, upib, marker) %>%
      mutate(upia = paste0(upia, "-A"), upib = paste0(upib, "-B")) %>%
      as_tbl_graph(directed = FALSE)

    # Add node type as node attribute
    g <- g %>%
      # Replace node attributes with name / node_type
      mutate(!!! g %>%
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
      new_edges <- g %E>% as_tibble() %>% select(from, to) %>% distinct()

      # Replace with collapsed edges and return a simple undirected graph
      g <- g %E>%
        select(-marker) %>%
        filter(FALSE) %>% # Remove all edges
        bind_edges(new_edges) %N>% # Add distinct edges and reactivate nodes
        select(name, node_type)
    }

    # Return results
    attr(g, "type") <- "bipartite"
    if (add_markers) {
      return(list(graph = g, counts = node_cntMatrix, counts_edges = edge_cntMatrix))
    } else {
      return(g)
    }
  })

  return(edge_table_split[cell_ids])

}


#' Load linegraph from edgelist
#'
#' @param arrow_data An arrow dataset
#' @param cell_id A character vector of cell IDs
#'
#' @import dplyr
#' @importFrom tidygraph as_tbl_graph
#'
#' @noRd
.load_as_linegraph <- function (
  arrow_data,
  cell_ids,
  add_markers = TRUE
) {

  # Start by loading graph as bipartite
  g_list <- .load_as_bipartite(arrow_data, cell_ids = cell_ids, add_markers = add_markers)

  g_list <- lapply(g_list, function(g) {

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
    if (add_markers) {
      # For the linegraph, the node counts are the same as
      # the edge counts of the bipartite graph
      return(list(graph = g_line, counts = g$counts_edges))
    } else {
      return(g_line)
    }
  })

  return(g_list)

}


#' Load A-node projected graph from edgelist
#'
#' @param arrow_data An arrow dataset
#' @param cell_id A cell ID
#'
#' @import dplyr
#' @importFrom tidygraph as_tbl_graph
#'
#' @noRd
.load_as_anode <- function (
  arrow_data,
  cell_ids,
  add_markers = TRUE
) {

  # Fetch edgelist from parquet file and group by component
  edge_table <- arrow_data %>%
    filter(component %in% cell_ids) %>%
    collect() %>%
    group_by(component)

  # Split into list with one component per element
  edge_table_split <- edge_table %>%
    group_split() %>%
    setNames(nm = edge_table %>% dplyr::group_data() %>% pull(component))

  edge_table_split <- lapply(edge_table_split, function(edge_table) {

    # Anode projection
    g_anode <- edge_table %>%
      select(upia, upib, marker) %>%
      left_join(edge_table,
                by = c("upib"),
                relationship = "many-to-many",
                suffix = c("1", "2")) %>%
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
        tidyr::pivot_wider(id_cols = upia, names_from = marker, values_from = n, values_fill = 0)

      cntMatrix <- as(markers_wide[, 2:ncol(markers_wide)] %>% as.matrix(), "dgCMatrix")

      # Sort matrix rows to match graph names
      nds <- g_anode %>%
        pull(name)
      cntMatrix <- cntMatrix[match(nds, markers_wide$upia), ]
      rownames(cntMatrix) <- nds
    }

    # Return results
    attr(g_anode, "type") <- "Anode"
    if (add_markers) {
      return(list(graph = g_anode, counts = cntMatrix))
    } else {
      return(g)
    }
  })

  # Return list with the correct order of components
  return(edge_table_split[cell_ids])

}

