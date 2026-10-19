#' Calculate Sequencing Saturation
#'
#' This function calculates the sequencing saturation of a sample or a graph component,
#' which can be applied to unique reads, nodes, or edges.
#'
#' @details
#' The sequencing saturation is calculated using the formula:
#' \deqn{S = 100 \times \left(1 - \frac{E}{R}\right)}
#' where:
#' - \eqn{S} is the saturation (as a percentage),
#' - \eqn{E} is the number of unique reads or graph elements (e.g., nodes, edges), and
#' - \eqn{R} is the number of reads or supporting elements in the graph.
#'
#' This function can be used for calculating the saturation of:
#' - Reads: Number of deduplicated reads vs. reads
#' - Nodes: Number of nodes vs. reads
#' - Edges: Number of edges vs. reads
#' - Other graph elements: Adjust the input accordingly
#'
#' @param graph_elements The number of graph elements (e.g., nodes, edges).
#' @param graph_reads The number of reads or supporting elements.
#'
#' @rdname sequencing_saturation
#' @importFrom tidygraph group_components
#'
#' @return The sequencing saturation of the graph expressed as a percentage.
#'
#' @examples
#'
#' # For a graph with 300 unique reads, 100 nodes, and 200 edges,
#' # sequenced at 400 total reads
#'
#' # Read sequencing saturation
#' sequencing_saturation(300, 400)
#'
#' # Node sequencing saturation
#' sequencing_saturation(100, 400)
#'
#' # Edge sequencing saturation
#' sequencing_saturation(200, 400)
#'
#' @export
#'
sequencing_saturation <- function(
  graph_elements,
  graph_reads
) {
  if (any(graph_elements > graph_reads)) {
    cli_warn("The number of graph elements should not exceed the number of reads.")
  }
  assert_vector(graph_elements, "numeric", n = 1)
  assert_vector(graph_reads, "numeric", n = 1)

  return(100 * (1 - (graph_elements / graph_reads)))
}


#' Simulate Sequencing Saturation Curve
#'
#' This function simulates the effect of lower read depth on the sequencing
#' saturation of a PNA sample. This can be used to create a saturation curve for
#' a given sample. The function iteratively downsamples the edgelist for a small
#' number of components and records the number of edges, proteins, and reads. The
#' saturation of the graph is calculated at each step.
#'
#' @param edgelist A tibble containing the edgelist with the following columns:
#' - `component`: The component number.
#' - `umi1`: The UMI of the first node
#' - `umi2`: The UMI of the second node.
#' - `read_count`: The number of reads supporting the edge.
#' @param sample_fracs A vector of sample fractions to downsample the edgelist.
#' The resulting sizes will not be exact, since some parts of the graph may be
#' disconnected upon downsampling. The function will keep the largest connected
#' component.
#' @param n_comps The number of components to sample from the edgelist.
#'
#' @rdname SequenceSaturationCurve
#'
#' @return A tibble with the following columns:
#' - `sample_size`: The number of reads in the downsampled edgelist.
#' - `sample_frac`: The fraction of reads in the downsampled edgelist.
#' - `graph_edges`: The number of edges in the downsampled graph.
#' - `graph_proteins`: The number of proteins in the downsampled graph.
#' - `graph_reads`: The number of reads in the downsampled graph.
#' - `graph_node_saturation`: The sequencing saturation of the graph based on the
#' number of proteins.
#' - `graph_edge_saturation`: The sequencing saturation of the graph based on the
#' number of edges.
#'
#' @examples
#'
#' library(dplyr)
#' library(ggplot2)
#' # Here we are reformatting an MPX edgelist to match the expected input
#' # for a PNA edgelist
#'
#' # Load the edgelist, and rename to match PNA input
#' edgelist <-
#'   ReadMPX_edgelist(system.file("extdata/five_cells", "five_cells.pxl",
#'     package = "pixelatorR"
#'   )) %>%
#'   rename(umi1 = upia, umi2 = upib, read_count = count)
#'
#' set.seed(37)
#' seqsat <- SequenceSaturationCurve(edgelist,
#'   sample_fracs = c(1, 0.75, 0.5, 0.25),
#'   n_comps = 2L
#' )
#'
#' # Calculate the mean sequencing saturation for each sample fraction
#' seqsat_mean <-
#'   seqsat %>%
#'   group_by(sample_frac) %>%
#'   summarise(
#'     mean_node_saturation = mean(graph_node_saturation),
#'     mean_edge_saturation = mean(graph_edge_saturation)
#'   )
#'
#' ggplot(seqsat_mean, aes(x = sample_frac, y = mean_node_saturation)) +
#'   geom_line() +
#'   labs(
#'     title = "Sequencing Saturation Curve",
#'     x = "Sample Fraction",
#'     y = "Node Saturation (%)"
#'   ) +
#'   theme_minimal()
#'
#' ggplot(seqsat_mean, aes(x = sample_frac, y = mean_edge_saturation)) +
#'   geom_line() +
#'   labs(
#'     title = "Sequencing Saturation Curve",
#'     x = "Sample Fraction",
#'     y = "Node Saturation (%)"
#'   ) +
#'   theme_minimal()
#'
#' @export
#'
SequenceSaturationCurve <- function(
  edgelist,
  sample_fracs = rev(seq(0.1, 1, 0.1)),
  n_comps = 10L
) {
  # This function simulates the effect of lower read depth on the
  # sequencing saturation of a sample. This can be used to create
  # a saturation curve for a given sample.

  assert_class(edgelist, "data.frame")
  assert_vector(sample_fracs, "numeric")
  assert_single_value(n_comps, type = "integer")
  if (!all(order(sample_fracs) == rev(seq_along(sample_fracs)))) {
    cli::cli_abort(
      c("x" = "{.var sample_fracs} must be in descending order")
    )
  }
  if (n_comps <= 0) {
    cli::cli_abort(
      c("x" = "{.var n_comps} must be greater than 0")
    )
  }

  # Read edgelist
  el_tot <-
    edgelist %>%
    filter(component %in% sample(unique(component), n_comps)) %>%
    select(component, umi1, umi2, read_count)

  tot_res <-
    el_tot %>%
    group_by(component) %>%
    do({
      el_init <-
        mutate(., edge = row_number()) %>%
        select(-component)

      index_init <-
        el_init %>%
        select(read_count) %>%
        pull(1) %>%
        as.numeric() %>%
        rep(seq_along(.), times = .)

      sample_sizes <-
        round(sample_fracs * length(index_init))

      full_graph_proteins <- n_distinct(c(el_init$umi1, el_init$umi2))

      res <- tibble()

      # Iteratively downsample the edgelist and record the number of
      # edges, proteins, and reads.
      for (i in seq_along(sample_sizes)) {
        if (i == 1) {
          index <- index_init
          el <- el_init
        } else {
          index <- sample(index_init, size = sample_sizes[i])
          el <- el_init[unique(index), ]

          # Create graph and remove all but the largest component
          el_graph <-
            el %>%
            as_tbl_graph() %>%
            mutate(component = group_components()) %>%
            filter(component == 1)

          # Update index
          index_keep <-
            el_graph %E>%
            pull(edge)

          index <- index[index %in% index_keep]
          el <- el[match(unique(index), el$edge), ]
        }
        graph_proteins <- n_distinct(c(el$umi1, el$umi2))
        res <-
          res %>%
          bind_rows(
            tibble(
              sample_size = sample_sizes[i],
              sample_frac = sample_fracs[i],
              graph_edges = nrow(el),
              graph_proteins = graph_proteins,
              graph_reads = length(index),
              graph_stability = graph_proteins / full_graph_proteins
            )
          )
      }

      res <-
        res %>%
        mutate(
          graph_node_saturation = sequencing_saturation(graph_proteins, graph_reads),
          graph_edge_saturation = sequencing_saturation(graph_edges, graph_reads)
        )

      res
    }) %>%
    ungroup()

  return(tot_res)
}

#' Compute approximate edge saturation
#'
#' Computes an approximate edge saturation for each component in the edgelist
#' using SQL queries on a DuckDB connection. Edge saturation measures how many
#' unique edges have been observed relative to the total number of reads.
#'
#' @details
#' Edge saturation is calculated as:
#' \deqn{S_{edge} = 1 - \frac{E}{R}}
#' where \eqn{E} is the number of unique edges and \eqn{R} is the total read count.
#'
#' Higher saturation values indicate that additional sequencing reads are unlikely
#' to reveal new edges, suggesting the sample has been sequenced to near-saturation.
#'
#' @param db A `PixelDB` object with an active DuckDB connection.
#' @param components Optional character vector of component names to filter the
#'   edgelist. If `NULL` (default), all components are included.
#' @param table_name Optional name for the computed remote table. If provided,
#'   the result is materialized as a table in the DuckDB connection.
#' @param min_degree Integer; minimum full-depth node degree required for both
#'   endpoints of an edge to be included. Default is `1L` (all edges included).
#'   Set to `2L` to exclude dangling (degree-1) leaves.
#'
#' @return A lazy `tbl` (DuckDB query) with columns:
#' \describe{
#'   \item{component}{The component (cell) identifier.}
#'   \item{edges}{The number of unique edges.}
#'   \item{read_count}{The total read count for the component.}
#'   \item{edge_saturation}{The edge saturation value (0 to 1).}
#' }
#'
#' @seealso
#' \code{\link{approximate_node_saturation}} for node-based saturation.
#' \code{\link{approximate_saturation_curve}} for computing saturation at
#' multiple downsampling fractions.
#'
#' @export
#'
approximate_edge_saturation <- function(
  db,
  components = NULL,
  table_name = NULL,
  min_degree = 1L
) {
  assert_class(db, "PixelDB")
  assert_vector(components, type = "character", allow_null = TRUE, n = 1)
  assert_single_value(table_name, type = "string", allow_null = TRUE)
  assert_single_value(min_degree, type = "integer")
  if (min_degree < 1L) {
    cli::cli_abort(
      c("x" = "{.var min_degree} must be >= 1, but got {.val {min_degree}}.")
    )
  }

  el_name <- .prepare_saturation_views_duckdb(
    db,
    components = components,
    min_degree = min_degree
  )

  edge_saturation <- tbl(db$.__enclos_env__$private$con, el_name) %>%
    .compute_saturation_duckdb(
      type = "edges",
      table_name = table_name
    )
  return(edge_saturation)
}

#' Compute approximate node saturation
#'
#' Computes an approximate node saturation for each component in the edgelist
#' using SQL queries on a DuckDB connection. Node saturation measures how many
#' unique nodes (UMIs) have been observed relative to the total number of reads.
#'
#' @details
#' Node saturation is calculated as:
#' \deqn{S_{node} = 1 - \frac{N}{R}}
#' where \eqn{N} is the number of unique nodes and \eqn{R} is the total read count
#' (counted twice per edge, once for each endpoint).
#'
#' This function creates a temporary view called `nodes` in the DuckDB connection
#' that aggregates node statistics (read count and degree) from the edgelist.
#'
#' Higher saturation values indicate that additional sequencing reads are unlikely
#' to reveal new nodes, suggesting the sample has been sequenced to near-saturation.
#'
#' @param db A `PixelDB` object with an active DuckDB connection.
#' @param components Optional character vector of component names to filter the
#'   edgelist. If `NULL` (default), all components are included.
#' @param table_name Optional name for the computed remote table. If provided,
#'   the result is materialized as a table in the DuckDB connection.
#' @param node_reads_multiplier Numeric; multiplier for the read count when calculating
#' node saturation. Default is `2`, reflecting that each edge contributes to two nodes.
#' @param min_degree Integer; minimum full-depth node degree required for nodes
#'   (and incident edges) to be included. Default is `1L` (all nodes included).
#'   Set to `2L` to exclude dangling (degree-1) leaves.
#'
#' @return A lazy `tbl` (DuckDB query) with columns:
#' \describe{
#'   \item{component}{The component (cell) identifier.}
#'   \item{nodes}{The number of unique nodes.}
#'   \item{read_count}{The total read count for the component scaled by
#'     `node_reads_multiplier` (i.e., twice the sequencing read count when
#'     `node_reads_multiplier = 2`, and equal to the sequencing read count when
#'     `node_reads_multiplier = 1`).}
#'   \item{node_saturation}{The node saturation value (0 to 1).}
#' }
#'
#' @seealso
#' \code{\link{approximate_edge_saturation}} for edge-based saturation.
#' \code{\link{approximate_saturation_curve}} for computing saturation at
#' multiple downsampling fractions.
#'
#' @export
#'
approximate_node_saturation <- function(
  db,
  components = NULL,
  table_name = NULL,
  node_reads_multiplier = 2L,
  min_degree = 1L
) {
  assert_class(db, "PixelDB")
  assert_vector(components, type = "character", allow_null = TRUE, n = 1)
  assert_single_value(table_name, type = "string", allow_null = TRUE)
  assert_single_value(node_reads_multiplier, type = "integer")
  if (!node_reads_multiplier %in% c(1, 2)) {
    cli::cli_abort(
      c("x" = "{.var node_reads_multiplier} must be either {.val 1}
               or {.val 2}, but got {.val {node_reads_multiplier}}.")
    )
  }
  assert_single_value(min_degree, type = "integer")
  if (min_degree < 1L) {
    cli::cli_abort(
      c("x" = "{.var min_degree} must be >= 1, but got {.val {min_degree}}.")
    )
  }

  .prepare_saturation_views_duckdb(
    db,
    components = components,
    min_degree = min_degree
  )

  # If node_reads_multiplier is 1, we treat the total component read count for
  # nodes as the same as for edges (half of the summed endpoint read counts).
  node_saturation <- tbl(db$.__enclos_env__$private$con, "nodes") %>%
    .compute_saturation_duckdb(
      type = "nodes",
      table_name = table_name,
      reads_multiplier = node_reads_multiplier / 2
    )

  return(node_saturation)
}

#' Prepare filtered edgelist and nodes views in DuckDB
#'
#' Filters the edgelist by `components` (if provided) and by `min_degree`
#' (retaining edges whose endpoints both have full-depth degree >= `min_degree`),
#' and creates the `nodes` temporary view from the resulting edgelist.
#'
#' @param db A `PixelDB` object with an active DuckDB connection.
#' @param components Optional character vector of component names to include.
#' @param min_degree Integer; minimum full-depth node degree to retain.
#' @param filter_edges Logical; if `TRUE` (default) and `min_degree > 1L`,
#'   filter the edgelist to edges whose endpoints both have full-depth
#'   `degree >= min_degree`. If `FALSE`, build `nodes` from the unfiltered
#'   edgelist and filter `nodes` by `degree >= min_degree` directly (avoiding
#'   self-joins on the edgelist when computing dynamic downsampling curves).
#'
#' @return Character string giving the name of the active edgelist table/view.
#'
#' @noRd
.prepare_saturation_views_duckdb <- function(
  db,
  components = NULL,
  min_degree = 1L,
  filter_edges = TRUE
) {
  db$check_connection()
  con <- db$.__enclos_env__$private$con

  if (is.null(components)) {
    el_name <- "edgelist"
  } else {
    if (!all(components %in% (db$counts() %>% colnames()))) {
      cli::cli_abort(
        c("x" = "All {.var components} must be present in the PXL file.")
      )
    }
    el_name <- "edgelist_modified"
    .filter_edgelist_duckdb(db, components)
  }

  if (min_degree > 1L && filter_edges) {
    DBI::dbExecute(
      con,
      glue::glue(
        "
        CREATE OR REPLACE TEMPORARY VIEW edgelist_deg_filtered AS (
          WITH u1 AS (
            SELECT component, umi1
            FROM {el_name}
            GROUP BY component, umi1
            HAVING COUNT(*) >= {as.integer(min_degree)}
          ),
          u2 AS (
            SELECT component, umi2
            FROM {el_name}
            GROUP BY component, umi2
            HAVING COUNT(*) >= {as.integer(min_degree)}
          )
          SELECT e.component, e.umi1, e.umi2, e.read_count
          FROM {el_name} e
          INNER JOIN u1 ON e.component = u1.component AND e.umi1 = u1.umi1
          INNER JOIN u2 ON e.component = u2.component AND e.umi2 = u2.umi2
        )
        "
      )
    )
    el_name <- "edgelist_deg_filtered"
  }

  having_clause <- if (min_degree > 1L && !filter_edges) {
    glue::glue("HAVING COUNT(*) >= {as.integer(min_degree)}")
  } else {
    ""
  }

  # Summarize node stats
  # Note that by default, each read supports two nodes (one at each endpoint),
  # so the total read count for nodes is 2x the total read count for edges.
  DBI::dbExecute(
    con,
    glue::glue(
      "
      CREATE OR REPLACE TEMPORARY VIEW nodes AS (
        SELECT component, CAST(SUM(read_count) AS INTEGER) AS read_count, CAST(COUNT(*) AS INTEGER) AS degree
        FROM {el_name}
        GROUP BY component, umi1
        {having_clause}
        UNION ALL
        SELECT component, CAST(SUM(read_count) AS INTEGER) AS read_count, CAST(COUNT(*) AS INTEGER) AS degree
        FROM {el_name}
        GROUP BY component, umi2
        {having_clause}
      )
      "
    )
  )

  return(el_name)
}

#' Filter edgelist by components
#'
#' Creates a temporary table `edgelist_modified` containing only edges from
#' the specified components.
#'
#' @param db A `PixelDB` object with an active DuckDB connection.
#' @param components Character vector of component names to include.
#'
#' @return Called for side effects. Creates a temporary table in the DuckDB
#'   connection.
#'
#' @noRd
.filter_edgelist_duckdb <- function(
  db,
  components
) {
  db$check_connection()
  components_formatted <- glue::glue_sql("{components}", .con = db$.__enclos_env__$private$con) %>%
    paste(collapse = ", ")
  DBI::dbExecute(
    db$.__enclos_env__$private$con,
    glue::glue(
      "CREATE OR REPLACE TEMPORARY TABLE edgelist_modified AS (
          SELECT component, umi1, umi2, read_count FROM edgelist WHERE component IN
          ({components_formatted})
        )"
    )
  )
}

#' Compute edge/node saturation from a lazy DuckDB table
#'
#' Computes saturation as `1 - (count of elements) / (sum of read counts)`
#' for each component.
#'
#' @param df_lazy A lazy `tbl` representing the edgelist or node table.
#' @param type Character; either `"edges"` or `"nodes"`.
#' @param table_name Optional name for materializing the result as a table.
#' @param reads_multiplier Numeric factor to scale the summed `read_count` before
#'   computing saturation (e.g., `0.5` when `node_reads_multiplier = 1` on a
#'   node table where each read appears at two endpoints).
#'
#' @return A lazy `tbl` with component, count, read_count, and saturation columns.
#'
#' @noRd
.compute_saturation_duckdb <- function(
  df_lazy,
  type = c("edges", "nodes"),
  table_name,
  reads_multiplier = 1
) {
  type <- match.arg(type, c("edges", "nodes"))
  type_singular <- stringr::str_replace(type, "s$", "")
  df_lazy %>%
    group_by(component) %>%
    summarize(
      !!sym(type) := as.integer(n()),
      read_count = round(sum(read_count, na.rm = TRUE) * reads_multiplier)
    ) %>%
    mutate(!!sym(paste0(type_singular, "_saturation")) := 1 - (!!sym(type) / read_count)) %>%
    compute(name = table_name, overwrite = TRUE)
}

#' Compute approximate saturation curve
#'
#' Computes an approximate saturation curve for nodes and edges in an edgelist
#' by estimating the expected number of retained elements at various downsampling
#' fractions. This provides insight into how sequencing depth affects saturation
#' without requiring actual downsampling.
#'
#' @details
#' ## Downsampling model
#'
#' The downsampling is performed probabilistically using expected values rather
#' than actual random sampling. For each edge with a given `read_count`, the
#' probability of the edge being retained at fraction `p` is:
#' \deqn{P(\text{edge retained}) = 1 - (1 - p)^{\text{read\_count}}}
#'
#' Similarly, a node is retained if at least one of its incident reads is
#' retained. With total incident reads \eqn{r_v = \sum_{e \ni v} r_e}, the
#' probability of node retention at fraction `p` is:
#' \deqn{P(\text{node retained}) = 1 - (1 - p)^{r_v}}
#'
#' ## Saturation calculation
#'
#' The saturation values are calculated using:
#' \deqn{S = 1 - \frac{N_{downsampled}}{R_{downsampled}}}
#'
#' where \eqn{N_{downsampled}} is the number of nodes or edges in the downsampled
#' graph, and \eqn{R_{downsampled}} is the number of reads supporting those elements.
#' Since each read supports two nodes (one at each endpoint),
#' \eqn{R_{downsampled, nodes} = 2 \times R_{downsampled, edges}}.
#'
#' ## Implementation notes
#'
#' This function calls \code{\link{approximate_node_saturation}} internally to
#' create the `nodes` temporary view in the DuckDB connection.
#'
#' @param db A `PixelDB` object with an active DuckDB connection.
#' @param fracs A numeric vector of fractions to downsample the edgelist.
#'   Values must be strictly between 0 and 1. Default is `seq(0.04, 0.96, by = 0.04)`.
#' @param components Optional character vector of component names to compute saturation for.
#' If `NULL` (default), all components are included.
#' @param node_reads_multiplier Numeric; multiplier for the read count when calculating
#' node saturation. Default is `2`, reflecting that each edge contributes to two nodes.
#' Optionally, set to `1` to use the same read count for nodes and edges, which may
#' be more appropriate in certain contexts.
#' @param min_degree Integer; minimum full-depth node degree required for nodes
#'   and incident edges to be included. Default is `1L` (all elements included).
#'   Set to `2L` to exclude dangling (degree-1) leaves.
#' @param verbose Logical; if `TRUE` (default), print progress messages.
#'
#' @return A tibble with the following columns:
#' \describe{
#'   \item{component}{The component (cell) identifier.}
#'   \item{read_count}{The estimated number of reads after downsampling.}
#'   \item{p}{The downsampling fraction.}
#'   \item{nodes}{The estimated number of nodes retained.}
#'   \item{edges}{The estimated number of edges retained.}
#'   \item{degree}{The estimated average node degree.}
#'   \item{node_saturation}{The node saturation at the given fraction.}
#'   \item{edge_saturation}{The edge saturation at the given fraction.}
#' }
#'
#' @seealso
#' \code{\link{approximate_node_saturation}} and \code{\link{approximate_edge_saturation}}
#' for computing saturation at full depth.
#' \code{\link{lcc_curve}} for computing largest connected component sizes at
#' different downsampling fractions.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' library(dplyr)
#'
#' pxl_file <- minimal_pna_pxl_file()
#' db <- PixelDB$new(pxl_file)
#'
#' sat_curve <- approximate_saturation_curve(db, fracs = seq(0.1, 0.9, by = 0.1))
#'
#' # Plot node saturation curve
#' sat_curve %>%
#'   tidyr::pivot_longer(
#'     cols = c(node_saturation, edge_saturation),
#'     names_to = "type",
#'     values_to = "saturation"
#'   ) %>%
#'   ggplot(aes(read_count, saturation, color = component)) +
#'   geom_line() +
#'   geom_point() +
#'   labs(x = "Reads", y = "Saturation") +
#'   theme_minimal() +
#'   facet_grid(~type) +
#'   scale_y_continuous(limits = c(0, 1))
#'
#' db$close()
#' }
#'
#' @export
#'
approximate_saturation_curve <- function(
  db,
  fracs = seq(0.04, 0.96, by = 0.04),
  components = NULL,
  node_reads_multiplier = 2L,
  min_degree = 1L,
  verbose = TRUE
) {
  assert_class(db, "PixelDB")
  assert_vector(fracs, type = "numeric", n = 1)
  assert_within_limits(fracs, limits = c(0, 1))
  assert_vector(components, type = "character", allow_null = TRUE, n = 1)
  assert_single_value(node_reads_multiplier, type = "integer")
  if (!node_reads_multiplier %in% c(1, 2)) {
    cli::cli_abort(
      c("x" = "{.var node_reads_multiplier} must be either {.val 1}
               or {.val 2}, but got {.val {node_reads_multiplier}}.")
    )
  }
  assert_single_value(min_degree, type = "integer")
  if (min_degree < 1L) {
    cli::cli_abort(
      c("x" = "{.var min_degree} must be >= 1, but got {.val {min_degree}}.")
    )
  }
  assert_single_value(verbose, type = "bool")

  # Prepare filtered edgelist and `nodes` temporary views in DuckDB
  el_name <- .prepare_saturation_views_duckdb(
    db,
    components = components,
    min_degree = min_degree,
    filter_edges = FALSE
  )

  if (verbose && check_global_verbosity()) {
    cli::cli_alert_info(
      "Computing downsampled nodes and edges for {.val {length(fracs)}} fractions..."
    )
  }

  # Pre-aggregate edges by (component, read_count) and nodes by
  # (component, read_count, degree) into compact histograms, then cross-join
  # with `fracs` via UNNEST. This reduces the number of POW() evaluations by
  # orders of magnitude compared to evaluating POW() per row per fraction.
  #
  # For min_degree == 1, a node or edge with `read_count` supporting reads is
  # retained at fraction `p` with probability `1 - (1 - p)^read_count`.
  # For min_degree >= 2, we dynamically estimate the probability that a node
  # has degree >= 2 at downsampling fraction `p` (i.e., at least 2 incident
  # edges retained at fraction `p`).
  fracs_sql <- paste(fracs, collapse = ", ")
  node_retention_sql <- if (min_degree == 1L) {
    "1.0 - POW(1.0 - f.p, n.read_count)"
  } else {
    paste0(
      "1.0 - POW(1.0 - f.p, n.read_count) - n.degree * (",
      "POW(1.0 - f.p, n.read_count * (n.degree - 1.0) / n.degree) - ",
      "POW(1.0 - f.p, n.read_count))"
    )
  }

  if (verbose && check_global_verbosity()) {
    cli::cli_alert_info(
      "Computing edge saturation, node saturation and average degree for downsampled graphs..."
    )
  }

  downsampled_saturation <- DBI::dbGetQuery(
    db$.__enclos_env__$private$con,
    glue::glue(
      "
      WITH fracs_tbl AS (
        SELECT CAST(UNNEST([{fracs_sql}]) AS DOUBLE) AS p
      ),
      edge_hist AS (
        SELECT component, read_count, COUNT(*) AS n_edges
        FROM {el_name}
        GROUP BY component, read_count
      ),
      node_hist AS (
        SELECT component, read_count, degree, COUNT(*) AS n_nodes
        FROM nodes
        GROUP BY component, read_count, degree
      ),
      edges_down AS (
        SELECT
          e.component,
          f.p,
          ROUND(SUM(e.n_edges * e.read_count) * f.p) AS read_count,
          SUM(e.n_edges * (1.0 - POW(1.0 - f.p, e.read_count))) AS edges
        FROM edge_hist e
        CROSS JOIN fracs_tbl f
        GROUP BY e.component, f.p
      ),
      nodes_down AS (
        SELECT
          n.component,
          f.p,
          SUM(n.n_nodes * ({node_retention_sql})) AS nodes
        FROM node_hist n
        CROSS JOIN fracs_tbl f
        GROUP BY n.component, f.p
      )
      SELECT
        ed.component,
        ed.read_count,
        ed.p,
        nd.nodes,
        ed.edges,
        2.0 * ed.edges / nd.nodes AS degree
      FROM edges_down ed
      JOIN nodes_down nd ON ed.component = nd.component AND ed.p = nd.p
      ORDER BY ed.component, ed.p
      "
    )
  ) %>%
    as_tibble() %>%
    mutate(
      node_saturation = 1 - (nodes / (read_count * node_reads_multiplier)),
      edge_saturation = 1 - (edges / read_count)
    ) %>%
    mutate(
      node_saturation = if_else(node_saturation < 0, NA_real_, node_saturation),
      edge_saturation = if_else(edge_saturation < 0, NA_real_, edge_saturation)
    )

  if (verbose && check_global_verbosity()) {
    cli::cli_alert_success("Finished!")
  }

  return(downsampled_saturation)
}

#' Compute marginal discovery rate and saturation
#'
#' Computes the marginal discovery rate (the instantaneous slope of the
#' discovery curve, i.e., expected new nodes and edges discovered per
#' additional sequencing read) and the corresponding marginal (Good-Turing)
#' saturation at full sequencing depth for each component in a `PixelDB`
#' object.
#'
#' @details
#' ## Marginal discovery rate (slope)
#'
#' Let \eqn{\mathbb{E}[U(p)]} be the expected number of unique elements (nodes
#' or edges) observed at downsampling fraction \eqn{p \in (0, 1]}, and let
#' \eqn{R(p) = p R} be the downsampled read count. The instantaneous slope
#' of the discovery curve with respect to sequencing reads is:
#' \deqn{\frac{d\,\mathbb{E}[U(p)]}{dR(p)} = \frac{1}{R} \frac{d\,\mathbb{E}[U(p)]}{dp}}
#'
#' Evaluated at full sequencing depth (\eqn{p = 1}):
#' - **Edges:** An edge with \eqn{r_e} reads has retention probability
#'   \eqn{1 - (1 - p)^{r_e}}, whose derivative at \eqn{p = 1} is \eqn{1} if
#'   \eqn{r_e = 1} and \eqn{0} if \eqn{r_e \ge 2}. Thus, the edge discovery
#'   rate at full depth is the fraction of reads supporting singleton edges
#'   (\eqn{E_1}):
#'   \deqn{\text{edge\_discovery\_rate} = \left.\frac{d\,\mathbb{E}[E]}{dR}\right|_{p=1} = \frac{E_1}{R}}
#' - **Nodes (`min_degree = 1`):** A node with total incident reads \eqn{r_v}
#'   has retention probability \eqn{1 - (1 - p)^{r_v}}, so at \eqn{p = 1}:
#'   \deqn{\text{node\_discovery\_rate} = \left.\frac{d\,\mathbb{E}[N]}{dR}\right|_{p=1} = \frac{N_1}{R}}
#'   where \eqn{N_1} is the number of singleton nodes (\eqn{r_v = 1}).
#' - **Core nodes (`min_degree = 2`):** A node is in the dynamic 2-core at
#'   fraction \eqn{p} if at least two of its incident edges are retained:
#'   \deqn{P(\text{deg}_v(p) \ge 2) = 1 + (d_v - 1)(1 - p)^{r_v} - \sum_{e \ni v} (1 - p)^{r_v - r_e}}
#'   Taking the derivative with respect to \eqn{p} at \eqn{p = 1}, the only
#'   terms with \eqn{r_v - r_e = 1} occur when a degree-2 node (\eqn{d_v = 2})
#'   has an incident singleton edge (\eqn{r_e = 1}), so:
#'   \deqn{\text{node\_discovery\_rate} =
#'     \left.\frac{d\,\mathbb{E}[N_{\text{deg} \ge 2}]}{dR}\right|_{p=1} =
#'     \frac{\sum_{v:\, d_v = 2} s_v}{R}}
#'   where \eqn{s_v \in \{0, 1, 2\}} is the number of incident singleton edges
#'   on degree-2 node \eqn{v}.
#'
#' Multiplying `node_discovery_rate` or `edge_discovery_rate` by `10000` gives
#' the expected number of new unique nodes or edges gained per 10,000 additional
#' sequencing reads.
#'
#' ## Marginal (Good-Turing) saturation
#'
#' The marginal saturation is the probability that the next observation lands
#' on an already-observed element rather than discovering a new one (since each
#' sequencing read observes 2 node endpoints, `node_discovery_rate` is divided
#' by `2`):
#' \deqn{S_{\text{marginal, edge}} = 1 - \text{edge\_discovery\_rate}}
#' \deqn{S_{\text{marginal, node}} = 1 - \frac{\text{node\_discovery\_rate}}{2}}
#'
#' @param db A `PixelDB` object with an active DuckDB connection.
#' @param components Optional character vector of component names to include.
#'   If `NULL` (default), all components are included.
#' @param table_name Optional name for materializing the result as a table in
#'   the DuckDB connection.
#' @param min_degree Integer (`1L` or `2L`); minimum node degree at the current
#'   depth for a node to be counted. Default is `1L` (all nodes). Set to `2L`
#'   to measure the discovery rate of dynamic degree >= 2 core nodes.
#'
#' @return A lazy `tbl` (DuckDB query) with columns:
#' \describe{
#'   \item{component}{The component (cell) identifier.}
#'   \item{read_count}{Total sequencing reads in the component.}
#'   \item{nodes}{Number of unique nodes (with degree >= `min_degree`).}
#'   \item{edges}{Number of unique edges in the component.}
#'   \item{node_discovery_rate}{Instantaneous slope \eqn{dN / dR} (new nodes per additional read).}
#'   \item{edge_discovery_rate}{Instantaneous slope \eqn{dE / dR} (new edges per additional read).}
#'   \item{marginal_node_saturation}{Marginal node saturation (\eqn{1 - (dN / dR) / 2}).}
#'   \item{marginal_edge_saturation}{Marginal (Good-Turing) edge saturation (\eqn{1 - dE / dR}).}
#' }
#'
#' @seealso
#' \code{\link{marginal_saturation_curve}} for computing marginal discovery
#' rates and saturation across downsampling fractions.
#' \code{\link{approximate_node_saturation}} and
#' \code{\link{approximate_edge_saturation}} for cumulative saturation.
#'
#' @export
#'
marginal_saturation <- function(
  db,
  components = NULL,
  table_name = NULL,
  min_degree = 1L
) {
  assert_class(db, "PixelDB")
  assert_vector(components, type = "character", allow_null = TRUE, n = 1)
  assert_single_value(table_name, type = "string", allow_null = TRUE)
  assert_single_value(min_degree, type = "integer")
  if (!min_degree %in% c(1L, 2L)) {
    cli::cli_abort(
      c("x" = "{.var min_degree} must be either {.val 1} or {.val 2},
               but got {.val {min_degree}}.")
    )
  }

  db$check_connection()
  con <- db$.__enclos_env__$private$con

  if (is.null(components)) {
    el_name <- "edgelist"
  } else {
    if (!all(components %in% (db$counts() %>% colnames()))) {
      cli::cli_abort(
        c("x" = "All {.var components} must be present in the PXL file.")
      )
    }
    el_name <- "edgelist_modified"
    .filter_edgelist_duckdb(db, components)
  }

  node_filter_sql <- if (min_degree == 1L) {
    "degree >= 1"
  } else {
    "degree >= 2"
  }

  node_singletons_sql <- if (min_degree == 1L) {
    "CASE WHEN read_count = 1 THEN 1 ELSE 0 END"
  } else {
    "CASE WHEN degree = 2 THEN singletons ELSE 0 END"
  }

  DBI::dbExecute(
    con,
    glue::glue(
      "
      CREATE OR REPLACE TEMPORARY VIEW marginal_saturation_view AS (
        WITH nodes_all AS (
          SELECT
            component,
            CAST(SUM(read_count) AS INTEGER) AS read_count,
            CAST(COUNT(*) AS INTEGER) AS degree,
            CAST(SUM(CASE WHEN read_count = 1 THEN 1 ELSE 0 END) AS INTEGER) AS singletons
          FROM {el_name}
          GROUP BY component, umi1
          UNION ALL
          SELECT
            component,
            CAST(SUM(read_count) AS INTEGER) AS read_count,
            CAST(COUNT(*) AS INTEGER) AS degree,
            CAST(SUM(CASE WHEN read_count = 1 THEN 1 ELSE 0 END) AS INTEGER) AS singletons
          FROM {el_name}
          GROUP BY component, umi2
        ),
        edge_stats AS (
          SELECT
            component,
            CAST(COUNT(*) AS INTEGER) AS edges,
            CAST(SUM(read_count) AS DOUBLE) AS read_count,
            CAST(SUM(CASE WHEN read_count = 1 THEN 1 ELSE 0 END) AS DOUBLE) AS edge_singletons
          FROM {el_name}
          GROUP BY component
        ),
        node_stats AS (
          SELECT
            component,
            CAST(SUM(CASE WHEN {node_filter_sql} THEN 1 ELSE 0 END) AS INTEGER) AS nodes,
            CAST(SUM({node_singletons_sql}) AS DOUBLE) AS node_marginal_num
          FROM nodes_all
          GROUP BY component
        )
        SELECT
          e.component,
          e.read_count,
          n.nodes,
          e.edges,
          n.node_marginal_num / e.read_count AS node_discovery_rate,
          e.edge_singletons / e.read_count AS edge_discovery_rate,
          1.0 - (n.node_marginal_num / (e.read_count * 2.0)) AS marginal_node_saturation,
          1.0 - (e.edge_singletons / e.read_count) AS marginal_edge_saturation
        FROM edge_stats e
        JOIN node_stats n ON e.component = n.component
      )
      "
    )
  )

  res <- tbl(con, "marginal_saturation_view")
  if (!is.null(table_name)) {
    res <- res %>% compute(name = table_name, overwrite = TRUE)
  }
  return(res)
}

#' Compute marginal discovery rate and saturation curve
#'
#' Computes the instantaneous slope of the node and edge discovery curves
#' (\eqn{dN / dR} and \eqn{dE / dR}) and the corresponding marginal
#' (Good-Turing) saturation across a vector of downsampling fractions `fracs`.
#'
#' @details
#' ## Marginal discovery rate across downsampling fractions
#'
#' At any downsampling fraction \eqn{p \in (0, 1]}, the expected number of
#' retained edges is \eqn{\mathbb{E}[E(p)] = \sum_{e=1}^{E} (1 - (1 - p)^{r_e})}.
#' Differentiating with respect to the downsampled read count \eqn{R(p) = p R}
#' gives the exact marginal edge discovery rate (new edges per additional read):
#' \deqn{\text{edge\_discovery\_rate}(p) =
#'   \frac{d\,\mathbb{E}[E(p)]}{dR(p)} =
#'   \frac{1}{R} \sum_{e=1}^{E} r_e (1 - p)^{r_e - 1}}
#'
#' For nodes:
#' - When `min_degree = 1L`, \eqn{\mathbb{E}[N(p)] = \sum_{v=1}^{N} (1 - (1 - p)^{r_v})},
#'   so the exact marginal node discovery rate is:
#'   \deqn{\text{node\_discovery\_rate}(p) =
#'     \frac{d\,\mathbb{E}[N(p)]}{dR(p)} =
#'     \frac{1}{R} \sum_{v=1}^{N} r_v (1 - p)^{r_v - 1}}
#' - When `min_degree = 2L`, a node has dynamic degree \eqn{\ge 2} at fraction
#'   \eqn{p} with exact probability:
#'   \deqn{P(\text{deg}_v(p) \ge 2) = 1 + (d_v - 1)(1 - p)^{r_v} - \sum_{e \ni v} (1 - p)^{r_v - r_e}}
#'   Differentiating with respect to \eqn{R(p) = p R} yields the exact marginal
#'   discovery rate of dynamic degree \eqn{\ge 2} core nodes:
#'   \deqn{\text{node\_discovery\_rate}(p) =
#'     \frac{1}{R} \sum_{v:\, d_v \ge 2} \left[
#'     \sum_{e \ni v} (r_v - r_e)(1 - p)^{r_v - r_e - 1} -
#'     (d_v - 1) r_v (1 - p)^{r_v - 1} \right]}
#'
#' @param db A `PixelDB` object with an active DuckDB connection.
#' @param fracs A numeric vector of downsampling fractions in \eqn{(0, 1]}.
#'   Default is `seq(0.04, 1.0, by = 0.04)`.
#' @param components Optional character vector of component names to include.
#'   If `NULL` (default), all components are included.
#' @param min_degree Integer (`1L` or `2L`); minimum dynamic node degree at
#'   fraction `p` for a node to be counted. Default is `1L` (all nodes). Set to
#'   `2L` to track dynamic degree >= 2 core nodes.
#' @param verbose Logical; if `TRUE` (default), print progress messages.
#'
#' @return A tibble with columns:
#' \describe{
#'   \item{component}{The component (cell) identifier.}
#'   \item{read_count}{Estimated number of reads at fraction `p`.}
#'   \item{p}{The downsampling fraction.}
#'   \item{nodes}{Expected number of unique nodes (with dynamic degree >= `min_degree`) at fraction `p`.}
#'   \item{edges}{Expected number of unique edges at fraction `p`.}
#'   \item{degree}{Estimated average node degree (\eqn{2 \times \text{edges} / \text{nodes}}).}
#'   \item{node_discovery_rate}{Instantaneous slope \eqn{dN / dR} at fraction `p` (new nodes per additional read).}
#'   \item{edge_discovery_rate}{Instantaneous slope \eqn{dE / dR} at fraction `p` (new edges per additional read).}
#'   \item{marginal_node_saturation}{Marginal node saturation at fraction `p` (\eqn{1 - (dN / dR) / 2}).}
#'   \item{marginal_edge_saturation}{Marginal (Good-Turing) edge saturation at fraction `p`.}
#' }
#'
#' @seealso
#' \code{\link{marginal_saturation}} for computing full-depth marginal
#' discovery rates and saturation as a lazy DuckDB table.
#' \code{\link{approximate_saturation_curve}} for cumulative saturation curves.
#'
#' @export
#'
marginal_saturation_curve <- function(
  db,
  fracs = seq(0.04, 1.0, by = 0.04),
  components = NULL,
  min_degree = 1L,
  verbose = TRUE
) {
  assert_class(db, "PixelDB")
  assert_vector(fracs, type = "numeric", n = 1)
  if (any(fracs <= 0 | fracs > 1)) {
    cli::cli_abort(
      c("x" = "All values in {.var fracs} must be in (0, 1].")
    )
  }
  assert_vector(components, type = "character", allow_null = TRUE, n = 1)
  assert_single_value(min_degree, type = "integer")
  if (!min_degree %in% c(1L, 2L)) {
    cli::cli_abort(
      c("x" = "{.var min_degree} must be either {.val 1} or {.val 2},
               but got {.val {min_degree}}.")
    )
  }
  assert_single_value(verbose, type = "bool")

  db$check_connection()
  con <- db$.__enclos_env__$private$con

  if (is.null(components)) {
    el_name <- "edgelist"
  } else {
    if (!all(components %in% (db$counts() %>% colnames()))) {
      cli::cli_abort(
        c("x" = "All {.var components} must be present in the PXL file.")
      )
    }
    el_name <- "edgelist_modified"
    .filter_edgelist_duckdb(db, components)
  }

  if (verbose && check_global_verbosity()) {
    cli::cli_alert_info(
      "Computing marginal discovery rates and saturation for {.val {length(fracs)}} fractions..."
    )
  }

  fracs_sql <- paste(fracs, collapse = ", ")

  if (min_degree == 1L) {
    query_sql <- glue::glue(
      "
      WITH fracs_tbl AS (
        SELECT CAST(UNNEST([{fracs_sql}]) AS DOUBLE) AS p
      ),
      nodes_all AS (
        SELECT component, CAST(SUM(read_count) AS INTEGER) AS read_count
        FROM {el_name}
        GROUP BY component, umi1
        UNION ALL
        SELECT component, CAST(SUM(read_count) AS INTEGER) AS read_count
        FROM {el_name}
        GROUP BY component, umi2
      ),
      edge_hist AS (
        SELECT component, read_count, COUNT(*) AS n_edges
        FROM {el_name}
        GROUP BY component, read_count
      ),
      node_hist AS (
        SELECT component, read_count, COUNT(*) AS n_nodes
        FROM nodes_all
        GROUP BY component, read_count
      ),
      edges_down AS (
        SELECT
          e.component,
          f.p,
          CAST(SUM(e.n_edges * e.read_count) AS DOUBLE) AS full_reads,
          ROUND(SUM(e.n_edges * e.read_count) * f.p) AS read_count,
          SUM(e.n_edges * (1.0 - POW(1.0 - f.p, e.read_count))) AS edges,
          SUM(e.n_edges * e.read_count * POW(1.0 - f.p, e.read_count - 1)) AS dE_dp
        FROM edge_hist e
        CROSS JOIN fracs_tbl f
        GROUP BY e.component, f.p
      ),
      nodes_down AS (
        SELECT
          n.component,
          f.p,
          SUM(n.n_nodes * (1.0 - POW(1.0 - f.p, n.read_count))) AS nodes,
          SUM(n.n_nodes * n.read_count * POW(1.0 - f.p, n.read_count - 1)) AS dN_dp
        FROM node_hist n
        CROSS JOIN fracs_tbl f
        GROUP BY n.component, f.p
      )
      SELECT
        ed.component,
        ed.read_count,
        ed.p,
        nd.nodes,
        ed.edges,
        2.0 * ed.edges / nd.nodes AS degree,
        nd.dN_dp / ed.full_reads AS node_discovery_rate,
        ed.dE_dp / ed.full_reads AS edge_discovery_rate
      FROM edges_down ed
      JOIN nodes_down nd ON ed.component = nd.component AND ed.p = nd.p
      ORDER BY ed.component, ed.p
      "
    )
  } else {
    query_sql <- glue::glue(
      "
      WITH fracs_tbl AS (
        SELECT CAST(UNNEST([{fracs_sql}]) AS DOUBLE) AS p
      ),
      u1 AS (
        SELECT component, umi1, CAST(SUM(read_count) AS INTEGER) AS read_count, CAST(COUNT(*) AS INTEGER) AS degree
        FROM {el_name}
        GROUP BY component, umi1
        HAVING COUNT(*) >= 2
      ),
      u2 AS (
        SELECT component, umi2, CAST(SUM(read_count) AS INTEGER) AS read_count, CAST(COUNT(*) AS INTEGER) AS degree
        FROM {el_name}
        GROUP BY component, umi2
        HAVING COUNT(*) >= 2
      ),
      nodes_deg2 AS (
        SELECT component, read_count, degree FROM u1
        UNION ALL
        SELECT component, read_count, degree FROM u2
      ),
      edge_hist AS (
        SELECT component, read_count, COUNT(*) AS n_edges
        FROM {el_name}
        GROUP BY component, read_count
      ),
      node_hist_2 AS (
        SELECT component, read_count, degree, COUNT(*) AS n_nodes
        FROM nodes_deg2
        GROUP BY component, read_count, degree
      ),
      r_other_hist AS (
        SELECT component, r_other, COUNT(*) AS n_inc
        FROM (
          SELECT e.component, u1.read_count - e.read_count AS r_other
          FROM {el_name} e
          JOIN u1 ON e.component = u1.component AND e.umi1 = u1.umi1
          UNION ALL
          SELECT e.component, u2.read_count - e.read_count AS r_other
          FROM {el_name} e
          JOIN u2 ON e.component = u2.component AND e.umi2 = u2.umi2
        )
        GROUP BY component, r_other
      ),
      edges_down AS (
        SELECT
          e.component,
          f.p,
          CAST(SUM(e.n_edges * e.read_count) AS DOUBLE) AS full_reads,
          ROUND(SUM(e.n_edges * e.read_count) * f.p) AS read_count,
          SUM(e.n_edges * (1.0 - POW(1.0 - f.p, e.read_count))) AS edges,
          SUM(e.n_edges * e.read_count * POW(1.0 - f.p, e.read_count - 1)) AS dE_dp
        FROM edge_hist e
        CROSS JOIN fracs_tbl f
        GROUP BY e.component, f.p
      ),
      nodes_down_base AS (
        SELECT
          n.component,
          f.p,
          SUM(n.n_nodes * (1.0 + (n.degree - 1.0) * POW(1.0 - f.p, n.read_count))) AS nodes_part1,
          SUM(-n.n_nodes * (n.degree - 1.0) * n.read_count * POW(1.0 - f.p, n.read_count - 1)) AS dN_dp_part1
        FROM node_hist_2 n
        CROSS JOIN fracs_tbl f
        GROUP BY n.component, f.p
      ),
      nodes_down_inc AS (
        SELECT
          r.component,
          f.p,
          SUM(r.n_inc * POW(1.0 - f.p, r.r_other)) AS nodes_part2,
          SUM(r.n_inc * r.r_other * POW(1.0 - f.p, r.r_other - 1)) AS dN_dp_part2
        FROM r_other_hist r
        CROSS JOIN fracs_tbl f
        GROUP BY r.component, f.p
      )
      SELECT
        ed.component,
        ed.read_count,
        ed.p,
        b.nodes_part1 - i.nodes_part2 AS nodes,
        ed.edges,
        2.0 * ed.edges / (b.nodes_part1 - i.nodes_part2) AS degree,
        (b.dN_dp_part1 + i.dN_dp_part2) / ed.full_reads AS node_discovery_rate,
        ed.dE_dp / ed.full_reads AS edge_discovery_rate
      FROM edges_down ed
      JOIN nodes_down_base b ON ed.component = b.component AND ed.p = b.p
      JOIN nodes_down_inc i ON ed.component = i.component AND ed.p = i.p
      ORDER BY ed.component, ed.p
      "
    )
  }

  curve_df <- DBI::dbGetQuery(con, query_sql) %>%
    as_tibble() %>%
    mutate(
      marginal_node_saturation = 1 - (node_discovery_rate / 2),
      marginal_edge_saturation = 1 - edge_discovery_rate
    ) %>%
    mutate(
      marginal_node_saturation = if_else(marginal_node_saturation < 0, NA_real_, marginal_node_saturation),
      marginal_edge_saturation = if_else(marginal_edge_saturation < 0, NA_real_, marginal_edge_saturation)
    )

  if (verbose && check_global_verbosity()) {
    cli::cli_alert_success("Finished!")
  }

  return(curve_df)
}


#' Downsample edgelist and export to parquet files
#'
#' Filters an edgelist by probabilistically downsampling read counts and exports
#' the results to parquet files. This simulates the effect of sequencing at lower
#' depth.
#'
#' @details
#' For each edge with a given `read_count`, the probability of the edge being
#' retained at fraction `p` is:
#' \deqn{P(\text{edge retained}) = 1 - (1 - p)^{\text{read\_count}}}
#'
#' This means edges with higher read counts are more likely to be retained,
#' which accurately models the effect of sequencing at lower depth.
#'
#' The filtered edgelists are written to parquet files in `outdir`, one file
#' per fraction. File names follow the pattern `edgelist_001.parquet`,
#' `edgelist_002.parquet`, etc., corresponding to the order of fractions.
#'
#' @param pxl_file Path to the PXL file.
#' @param outdir Directory to write the output parquet files. Created if it
#'   does not exist.
#' @param components A character vector of component names to include. All
#'   components must be present in the PXL file.
#' @param fracs A numeric vector of fractions to downsample by. Values must
#'   be strictly between 0 and 1. Default is `seq(0.1, 0.9, 0.1)`.
#'
#' @return A tibble with columns:
#' \describe{
#'   \item{fracs}{The downsampling fraction.}
#'   \item{pq_files}{Path to the corresponding parquet file.}
#' }
#'
#' @seealso
#' \code{\link{lcc_sizes}} for computing LCC sizes from the output parquet files.
#' \code{\link{lcc_curve}} for a high-level wrapper that combines downsampling
#' and LCC computation.
#'
#' @export
#'
downsample_to_parquet <- function(
  pxl_file,
  outdir,
  components,
  fracs = seq(0.1, 0.9, 0.1)
) {
  assert_pxl_file(pxl_file)
  assert_single_value(outdir, type = "string")
  if (!fs::dir_exists(outdir)) {
    fs::dir_create(outdir)
  }
  assert_vector(fracs, type = "numeric", n = 1)
  assert_within_limits(fracs, limits = c(0, 1))

  # Establish a connection to the PXL file
  db <- PixelDB$new(pxl_file)
  on.exit({
    db$close()
    rm(db)
  })
  if (!all(components %in% (db$counts() %>% colnames()))) {
    cli::cli_abort("All {.var components} must be present in the PXL file")
  }
  comps_formatted <- glue::glue_sql("{components}", .con = db$.__enclos_env__$private$con) %>% paste(collapse = ",")

  pq_paths <- file.path(outdir, glue::glue("edgelist_{sprintf('%03d', seq_along(fracs))}.parquet"))
  for (f in pq_paths) {
    if (fs::file_exists(f)) {
      fs::file_delete(f)
    }
  }

  # Iterate over selected fractions
  for (i in seq_along(fracs)) {
    p <- fracs[i]
    pq_path <- pq_paths[i]
    DBI::dbExecute(
      db$.__enclos_env__$private$con,
      glue::glue(
        "
      COPY (
        SELECT q01.*
        FROM (
          SELECT q01.*, RANDOM() < pr AS keep
          FROM (
            SELECT
              umi1,
              umi2,
              read_count,
              component,
              1.0 - POW((1.0 - {p}), read_count) AS pr -- Your probability calculation
            FROM edgelist
            WHERE (component IN ({comps_formatted}))
          ) q01
        ) q01
        WHERE (keep)
      )
      TO '{pq_path}' (FORMAT PARQUET);
      "
      )
    )
  }

  return(tibble(fracs = fracs, pq_files = pq_paths))
}

#' Compute LCC sizes from downsampled edgelists
#'
#' Calculates the sizes of the largest connected components (LCCs) from a set
#' of parquet files containing downsampled edgelists.
#'
#' @details
#' ## Algorithm
#'
#' For each parquet file (representing a downsampling fraction), this function:
#' 1. Loads the edgelist into a DuckDB connection
#' 2. Constructs a property graph using the `duckpgq` extension
#' 3. Computes weakly connected components
#' 4. Identifies the largest component for each cell
#'
#' ## Performance
#'
#' The `duckpgq` extension enables efficient graph operations without loading
#' the full edgelist into memory. Parallel processing across fractions is
#' supported via `mc_cores`.
#'
#' @param df A tibble with columns:
#' \describe{
#'   \item{pq_files}{Paths to parquet files containing downsampled edgelists.}
#'   \item{fracs}{The corresponding downsampling fractions.}
#' }
#' Typically generated with \code{\link{downsample_to_parquet}}.
#' @param mc_cores Number of cores to use for parallel processing. Default is 1
#'   (sequential processing). Values > 1 use `parallel::mclapply`.
#'
#' @return A tibble with columns:
#' \describe{
#'   \item{component}{The component (cell) identifier.}
#'   \item{frac}{The downsampling fraction.}
#'   \item{n_nodes}{The number of nodes in the largest connected component.}
#' }
#'
#' @seealso
#' \code{\link{downsample_to_parquet}} for generating the input parquet files.
#' \code{\link{lcc_curve}} for a high-level wrapper.
#'
#' @export
#'
lcc_sizes <- function(
  df,
  mc_cores = 1
) {
  assert_class(df, "tbl_df")
  for (f in df$pq_files) {
    assert_file_exists(f)
    assert_file_ext(f, "parquet")
  }
  assert_within_limits(df$fracs, limits = c(0, 1))
  assert_single_value(mc_cores, type = "integer")

  lcc <- parallel::mclapply(seq_len(nrow(df)), function(i) {
    # Create a new DuckDB connection and load the duckpgq extension
    con <- DBI::dbConnect(duckdb::duckdb(bigint = "integer64"), dbdir = ":memory:")
    on.exit(DBI::dbDisconnect(con), add = TRUE)

    DBI::dbExecute(con, "INSTALL duckpgq FROM community")
    DBI::dbExecute(con, "LOAD duckpgq")

    # Create table from the current parquet file
    pq_path <- df$pq_files[i]
    DBI::dbExecute(con, glue::glue("CREATE OR REPLACE TABLE el AS SELECT * FROM '{pq_path}'"))

    # Create a node table from unique UMIs
    DBI::dbExecute(
      con,
      "
        CREATE OR REPLACE TEMPORARY TABLE umi_nodes AS
        SELECT component, umi1 AS id FROM el
        UNION
        SELECT component, umi2 AS id FROM el
      "
    )

    # Construct property graph from the node table and edgelist
    DBI::dbExecute(
      con,
      "
        CREATE OR REPLACE PROPERTY GRAPH my_umi_graph
        VERTEX TABLES (
            umi_nodes PROPERTIES (id) LABEL UMI
        )
        EDGE TABLES (
            el SOURCE KEY (umi1) REFERENCES umi_nodes (id)
                     DESTINATION KEY (umi2) REFERENCES umi_nodes (id)
                     LABEL CONNECTS
        );
      "
    )

    # Calculate weakly connected components
    DBI::dbExecute(
      con,
      "
        CREATE OR REPLACE TEMPORARY TABLE comps AS
        SELECT *
        FROM weakly_connected_component(my_umi_graph, UMI, CONNECTS);
      "
    )

    # Get the size of the largest connected components (LCC)
    # for each cell (component) in the edgelist
    lcc_per_component <- tbl(con, "comps") %>%
      left_join(tbl(con, "umi_nodes"), by = c("id" = "id")) %>%
      group_by(componentId, component) %>%
      count(name = "n_nodes") %>%
      mutate(n_nodes = as.integer(n_nodes)) %>%
      group_by(component) %>%
      arrange(desc(n_nodes)) %>%
      collect() %>%
      slice_head(n = 1)

    return(lcc_per_component)
  }, mc.cores = mc_cores)

  # Format results into a single tibble
  lcc <- lapply(seq_along(lcc), function(i) {
    lcc[[i]] %>%
      mutate(frac = df$fracs[i]) %>%
      select(component, frac, n_nodes)
  }) %>%
    bind_rows() %>%
    ungroup()

  return(lcc)
}


#' Compute LCC curve for downsampled edgelists
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Computes the sizes of the largest connected components (LCC) for edgelists
#' downsampled at multiple fractions. This is useful for analyzing graph stability
#' and determining the critical sequencing depth.
#'
#' @details
#' ## Graph stability analysis
#'
#' The LCC (largest connected component) measures graph connectivity. When
#' downsampling sequencing reads, edges are removed probabilistically, which
#' eventually causes the graph to fragment into smaller disconnected components.
#'
#' By computing LCC sizes at multiple downsampling fractions, you can observe
#' how graph connectivity degrades with reduced sequencing depth. Typically,
#' there is a sharp transition point where the LCC collapses, indicating the
#' minimum sequencing depth required to maintain graph structure.
#'
#' ## Workflow
#'
#' This function is a high-level wrapper that:
#' 1. Calls \code{\link{downsample_to_parquet}} to create downsampled edgelists
#' 2. Calls \code{\link{lcc_sizes}} to compute LCC sizes for each fraction
#' 3. Optionally cleans up temporary files
#'
#' ## Performance
#'
#' Computations are scalable thanks to efficient SQL queries. LCCs are computed
#' using the `duckpgq` extension for DuckDB, enabling fast graph operations
#' without loading the full edgelist into memory.
#'
#' For large datasets, it is recommended to compute LCC sizes on a subset of
#' components to minimize memory usage and computation time.
#'
#' @param pxl_file Path to the input PXL file containing the edgelist.
#' @param components A character vector of component names to include. All
#'   components must be present in the PXL file.
#' @param fracs A numeric vector of downsampling fractions. Values must be
#'   strictly between 0 and 1. Default is `seq(0.1, 1, by = 0.1)`.
#' @param outdir Optional output directory for downsampled parquet files.
#'   If `NULL` (default), files are saved to a temporary directory and deleted
#'   after computation.
#' @param mc_cores Number of cores for parallel processing. Default is 1
#'   (sequential). Values > 1 use `parallel::mclapply`.
#' @param verbose Logical; if `TRUE` (default), print progress messages.
#'
#' @return A tibble with columns:
#' \describe{
#'   \item{component}{The component (cell) identifier.}
#'   \item{frac}{The downsampling fraction.}
#'   \item{n_nodes}{The number of nodes in the largest connected component.}
#'   \item{read_count}{The estimated number of reads supporting the LCC at the given fraction.}
#' }
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' library(dplyr)
#' pxl_file <- minimal_pna_pxl_file()
#' components <- ReadPNA_counts(pxl_file) %>% colnames()
#'
#' # Select fractions based on average read depth
#' db <- PixelDB$new(pxl_file)
#' avg_reads <- db$cell_meta()$reads %>% mean()
#' fracs <- (10^seq(4, log10(avg_reads), length.out = 20))[-20] / avg_reads
#'
#' lcc_df <- lcc_curve(pxl_file, components, fracs = fracs)
#'
#' # Plot LCC sizes by fraction
#' ggplot(lcc_df, aes(frac, n_nodes, color = component)) +
#'   geom_point() +
#'   geom_line() +
#'   theme_bw() +
#'   guides(color = "none")
#'
#' # Compute graph stability (LCC / max theoretical nodes)
#' nodesat <- approximate_node_saturation(db) %>% collect()
#' db$close()
#'
#' lcc_df <- lcc_df %>%
#'   left_join(nodesat, by = "component")
#'
#' ggplot(lcc_df, aes(frac, n_nodes / nodes, color = component)) +
#'   geom_point() +
#'   geom_line() +
#'   scale_y_continuous(limits = c(0, 1)) +
#'   labs(y = "Graph stability (LCC / total nodes)") +
#'   theme_bw() +
#'   guides(color = "none")
#' }
#'
#' @seealso
#' \code{\link{downsample_to_parquet}} for the downsampling step.
#' \code{\link{lcc_sizes}} for computing LCC sizes from parquet files.
#' \code{\link{approximate_saturation_curve}} for saturation analysis without
#' actual downsampling.
#'
#' @export
#'
lcc_curve <- function(
  pxl_file,
  components,
  fracs = seq(0.1, 1, by = 0.1),
  outdir = NULL,
  mc_cores = 1,
  verbose = TRUE
) {
  duckdb_v <- utils::packageVersion("duckdb")
  if (utils::compareVersion(as.character(duckdb_v), "1.5.0") > 0) {
    cli_alert_warning(
      "This function is only tested with duckdb <= 1.5.0, but you have {.pkg duckdb} version {.val {duckdb_v}}
      installed. This may result in compatibility issues with {.pkg duckpgq}."
    )
  }

  assert_single_value(pxl_file, type = "string")
  assert_pxl_file(pxl_file)
  assert_vector(components, type = "character", n = 1)
  assert_single_value(outdir, "string", allow_null = TRUE)
  assert_vector(fracs, type = "numeric", allow_null = FALSE, n = 1)
  assert_within_limits(fracs, limits = c(0, 1))
  assert_single_value(mc_cores, type = "integer")
  if (!is.null(outdir)) {
    if (!fs::dir_exists(outdir)) {
      cli::cli_abort("{.var outdir} {.val {outdir}} does not exist")
    }
    clean_up <- FALSE
  } else {
    outdir <- fs::file_temp() %>% stringr::str_replace("file", "dir")
    fs::dir_create(outdir)
    clean_up <- TRUE
  }

  if (verbose && check_global_verbosity()) {
    cli::cli_alert_info("Downsampling and exporting edgelists for {.val {length(components)}} components")
  }

  pq_files <- downsample_to_parquet(pxl_file, outdir, components, fracs)

  if (verbose && check_global_verbosity()) {
    cli::cli_alert_info("Calculating largest connected components (LCC) for each fraction")
  }

  lcc <- lcc_sizes(df = pq_files, mc_cores)

  # Add read counts to the LCC results for better interpretability
  db <- PixelDB$new(pxl_file)
  on.exit({
    db$close()
  })
  read_counts <- db$cell_meta() %>%
    select(read_count = reads_in_component) %>%
    rownames_to_column("component")
  lcc <- lcc %>%
    left_join(read_counts, by = "component") %>%
    mutate(read_count = as.integer(round(read_count * frac)))

  if (clean_up) {
    fs::dir_delete(outdir)
  }

  if (verbose && check_global_verbosity()) {
    cli::cli_alert_success("Finished!")
  }

  return(lcc)
}
