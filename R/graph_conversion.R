#' @include generics.R
NULL

#' @param cl A cluster object created by makeCluster, or an integer
#' to indicate number of child-processes (integer values are ignored
#' on Windows) for parallel evaluations. See Details on performance
#' in the documentation for \code{pbapply}. The default is NULL,
#' which means that no parallelization is used.
#'
#' @rdname graph-conversion
#' @method edgelist_to_simple_Anode_graph data.frame
#'
#' @examples
#'
#' library(pixelatorR)
#' library(tibble)
#'
#' pxl_file <- system.file("extdata/five_cells", "five_cells.pxl", package = "pixelatorR")
#'
#' # Load edgelist
#' el <- ReadMPX_arrow_edgelist(pxl_file)
#'
#' # Convert to tbl_df
#' el_tbl_df <- as_tibble(el)
#'
#' # data.frame method --------------------------------
#' # Create list of A-node projected tbl_graphs
#' anode_proj <- edgelist_to_simple_Anode_graph(el_tbl_df)
#'
#' @export
#'
edgelist_to_simple_Anode_graph.data.frame <- function(
  object,
  components = NULL,
  cl = NULL,
  verbose = TRUE,
  ...
) {
  # Check input parameters
  stopifnot(
    "edgelist must be a non-empty object" =
      (nrow(object) > 0),
    "One or several of 'upia', 'upib' are missing from edgelist" =
      all(c("upia", "upib") %in% colnames(object))
  )
  if (!is.null(components) && ("component" %in% colnames(object))) {
    stopifnot(
      "'components' must be a character vector" =
        is.character(components) && (length(components) > 0)
    )
    if (!all(components %in% (object %>% pull(component)))) {
      abort("Some 'component' IDs are missing from object")
    }
    # Filter components
    object <- object %>%
      filter(component %in% components)
  }

  # Split by cell/component
  if ("component" %in% colnames(object)) {
    edgelist_split <- object %>%
      group_by(component)
    components <- edgelist_split %>%
      group_data() %>%
      pull(component)
    edgelist_split <- edgelist_split %>%
      group_split()
  } else {
    edgelist_split <- list(object)
  }

  # Simplify edgelist
  if (verbose && check_global_verbosity()) {
    cli_alert_info("Simplifying edge list")
  }
  edgelist_split <- lapply(edgelist_split, function(edgelist) {
    edgelist <-
      edgelist %>%
      select(upia, upib) %>%
      distinct()
  })

  # Create Anode graph
  if (verbose && check_global_verbosity()) {
    cli_alert_info("Creating A-node projected graphs")
  }

  edgelist_split <- pblapply(edgelist_split, function(edgelist) {
    anode_graph <-
      edgelist %>%
      left_join(edgelist,
        by = c("upib"),
        relationship = "many-to-many",
        suffix = c("1", "2")
      ) %>%
      select(upia1, upia2) %>%
      distinct() %>%
      filter(upia1 < upia2)
    return(anode_graph)
  }, cl = cl)

  # Convert edge lists to tbl_graphs and add attribute
  anode_graphs <- lapply(edgelist_split, function(edgelist) {
    g <- as_tbl_graph(edgelist, directed = FALSE)
    attr(g, "type") <- "Anode"
    return(g)
  }) %>%
    set_names(nm = components)

  if (verbose && check_global_verbosity()) {
    cli_alert_success("Returning an A-node projected graphs")
  }
  return(anode_graphs)
}


#' @import rlang
#'
#' @rdname graph-conversion
#' @method edgelist_to_simple_Anode_graph FileSystemDataset
#'
#' @examples
#' # FileSystemDataset method -------------------------
#' # Create list of A-node projected tbl_graphs
#' anode_proj <- edgelist_to_simple_Anode_graph(el)
#'
#' @export
#'
edgelist_to_simple_Anode_graph.FileSystemDataset <- function(
  object,
  components = NULL,
  verbose = TRUE,
  ...
) {
  expect_duckdb()

  # Check input parameters
  stopifnot(
    "edgelist must be a non-empty object" =
      (nrow(object) > 0),
    "One or several of 'upia', 'upib' are missing from edgelist" =
      all(c("upia", "upib") %in% names(object))
  )
  if (!"component" %in% names(object)) {
    abort("Column 'component' is missing from edgelist")
  }

  object <- object %>% to_duckdb()

  if (!is.null(components)) {
    stopifnot(
      "'components' must be a character vector" =
        is.character(components) && (length(components) > 0)
    )
    if (!all(components %in% (object %>% pull(component)))) {
      abort("Some 'component' IDs are missing from object")
    }
    # Filter components
    object <- object %>%
      filter(component %in% components)
  }

  object <- object %>%
    select(upia, upib, component) %>%
    group_by(component)

  # Simplify edgelist
  if (verbose && check_global_verbosity()) {
    cli_alert_info("Simplifying edge list")
  }
  object <- object %>%
    distinct()

  # Add suffix
  object <- object %>%
    mutate(rn = paste0("_", row_number())) %>%
    mutate(upia = str_c(upia, rn)) %>%
    select(-rn)

  if (verbose && check_global_verbosity()) {
    cli_alert_info("Creating A-node projected graphs")
  }
  anode_graph <- object %>%
    left_join(
      y = object,
      by = c("upib", "component"),
      suffix = c("1", "2")
    ) %>%
    select(upia1, upia2, component) %>%
    group_by(component) %>%
    mutate(
      upia1 = str_sub(upia1, start = 1, end = 25),
      upia2 = str_sub(upia2, start = 1, end = 25)
    ) %>%
    distinct() %>%
    filter(upia1 < upia2) %>%
    collect()

  # Split into multiple graphs
  components <- anode_graph %>%
    group_data() %>%
    pull(component)
  anode_graphs <- anode_graph %>%
    group_split()

  # Convert to tbl_graph
  anode_graphs <- lapply(anode_graphs, function(anode_graph) {
    g <- as_tbl_graph(anode_graph %>% select(-component), directed = FALSE)
    attr(g, "type") <- "Anode"
    return(g)
  }) %>%
    set_names(nm = components)

  if (verbose && check_global_verbosity()) {
    cli_alert_success("Returning an A-node projected graphs")
  }
  return(anode_graphs)
}


#' Create a simple bipartite graph from an edgelist
#'
#' @param edgelist An object of class \code{tbl_graph}
#'
#' @return An object of class \code{tbl_graph} containing a bipartite graph
#'
#' @export
#'
edgelist_to_simple_bipart_graph <- function(
  edgelist
) {
  # Check component column
  if ("component" %in% names(edgelist)) {
    no_components <- length(unique(edgelist %>% pull(component)))
    if (no_components > 1) {
      abort(glue("Found {no_components} components in edgelist."))
    }
  }

  # Simplify edgelist
  edgelist <- edgelist %>%
    select(upia, upib, marker) %>%
    distinct()

  # Create bipartite, tidy graph
  bipart_graph <- edgelist %>%
    as_tbl_graph(directed = FALSE) %>%
    mutate(node_type = case_when(name %in% edgelist$upia ~ "A", TRUE ~ "B"))
  attr(bipart_graph, "type") <- "bipartite"
  return(bipart_graph)
}
