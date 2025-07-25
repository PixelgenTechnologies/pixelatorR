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
          # Each edge contributes to 2 nodes, hence each edge read contributes to 2
          # nodes. Therefore, the total reads supporting the nodes is double the
          # number of reads supporting edges.
          graph_node_saturation = sequencing_saturation(graph_proteins, graph_reads * 2),
          graph_edge_saturation = sequencing_saturation(graph_edges, graph_reads)
        )

      res
    }) %>%
    ungroup()

  return(tot_res)
}

#' Compute approximate edge saturation
#'
#' This function computes an approximate edge saturation for each component
#' in the edgelist. It estimates the theoretical maximum number of edges
#' using the Chao1 estimator and calculates the edge saturation as the ratio
#' of the actual number of edges to the theoretical maximum.
#'
#' @param db A `PixelDB` object.
#' @param table_name Optional name of the remote table.
#'
#' @return A `lazy_df` with the edge saturation.
#'
#' @export
#'
approximate_edge_saturation <- function(
  db,
  table_name = NULL
) {
  assert_class(db, "PixelDB")
  assert_single_value(table_name, type = "string", allow_null = TRUE)
  # Estimate chao1 edge saturation
  edge_saturation <- tbl(db$.__enclos_env__$private$con, "edgelist") %>%
    .compute_saturation_lazy(
      type = "edges",
      table_name = table_name
    )
  return(edge_saturation)
}

#' Compute approximate node saturation
#'
#' This function computes an approximate node saturation for each component
#' in the edgelist. It estimates the theoretical maximum number of nodes
#' using the Chao1 estimator and calculates the nodes saturation as the ratio
#' of the actual number of nodes to the theoretical maximum.
#'
#' @param db A `PixelDB` object.
#' @param table_name Optional name of the remote table.
#'
#' @return A `lazy_df` with the node saturation.
#'
#' @export
#'
approximate_node_saturation <- function(
  db,
  table_name = NULL
) {
  assert_class(db, "PixelDB")
  assert_single_value(table_name, type = "string", allow_null = TRUE)
  # Summarize node stats
  DBI::dbExecute(
    db$.__enclos_env__$private$con,
    "
    CREATE OR REPLACE TEMPORARY VIEW nodes AS (
      SELECT component, CAST(SUM(read_count) AS SMALLINT) AS read_count, CAST(COUNT(*) AS SMALLINT) AS degree
      FROM edgelist
      GROUP BY component, umi1
      UNION ALL
      SELECT component, CAST(SUM(read_count) AS SMALLINT) AS read_count, CAST(COUNT(*) AS SMALLINT) AS degree
      FROM edgelist
      GROUP BY component, umi2
    )
    "
  )

  node_saturation <- tbl(db$.__enclos_env__$private$con, "nodes") %>%
    .compute_saturation_lazy(
      type = "nodes",
      table_name = table_name
    )

  return(node_saturation)
}

#' Compute edge/node saturation using Chao1 estimator
#'
#' @noRd
#'
.compute_saturation_lazy <- function(
  df_lazy,
  type = c("edges", "nodes"),
  table_name
) {
  type <- match.arg(type, c("edges", "nodes"))
  type_singular <- stringr::str_replace(type, "s$", "")
  df_lazy %>%
  group_by(component) %>%
    summarize(
      !! sym(type) := as.integer(n()),
      f1 = sum(read_count == 1, na.rm = TRUE),
      f2 = sum(read_count == 2, na.rm = TRUE)
    ) %>%
    mutate(
      !! sym(paste0("theoretical_max_", type)) := !! sym(type) + (f1 * (f1 - 1)) / (2 * (f2 + 1)),
    ) %>%
    mutate(
      !! sym(paste0(type_singular, "_saturation")) :=
        !! sym(type) / !! sym(paste0("theoretical_max_", type))
    ) %>%
    select(all_of(c("component", type, paste0(type_singular, "_saturation"), paste0("theoretical_max_", type)))) %>%
    compute(name = table_name, overwrite = TRUE)
}

#' Compute approximate saturation curve
#'
#' This function computes an approximate saturation curve for nodes and edges
#' in an edgelist. The edgelist is downsampled to various fractions of the
#' total number of reads, and the number of remaining nodes and edges are calculated
#' for each fraction. The saturation values are then calculated relative to the
#' theoretical maximum number of nodes and edges using the Chao1 estimator. The
#' theoretical maximum is computed for each component in the full edgelist. Note
#' that Chao1 will give a lower bound for the theoretical maximum, hence the
#' saturation values are likely overestimated. The estimate will be more robust
#' if the sample was sequenced at a high depth.
#'
#' @param db A `PixelDB` object.
#' @param fracs A numeric vector of fractions to downsample the edgelist.
#' @param detailed If `TRUE`, the function will return ...
#' @param verbose Print messages
#'
#' @return A tibble with the saturation values
#'
#' @export
#'
approximate_saturation_curve <- function(
  db,
  fracs = seq(0.04, 0.96, by = 0.04),
  detailed = FALSE,
  verbose = TRUE
) {
  assert_class(db, "PixelDB")
  assert_vector(fracs, type = "numeric", n = 1)
  assert_within_limits(fracs, limits = c(0, 1))
  assert_single_value(detailed, type = "bool")

  if (verbose && check_global_verbosity()) {
    cli::cli_alert_info("Computing theoretical maxmimum for nodes and edges using the Chao1 estimator... ")
  }

  # Compute node and edge saturation per component
  node_saturation <- approximate_node_saturation(db, table_name = "node_saturation")
  edge_saturation <- approximate_edge_saturation(db, table_name = "edge_saturation")

  if (verbose && check_global_verbosity()) {
    cli::cli_alert_info("Computing downsampled nodes and edges for break points:\n {.val {fracs}}... ")
  }

  # Calculate expected fraction of edges remaining after downsampling
  mut_actual_p <- list()
  for (p in fracs) {
    col_name <- paste0("p_", p)
    # Probability of being removed: (1 - p) ^ read_count
    mut_actual_p[[col_name]] <- rlang::expr((1 - !!p) ^ read_count)
  }
  sum_est_edges <- list()
  for (col_name in paste0("p_", fracs)) {
    summary_col_name <- paste0(col_name, "_tot")
    # Total expected edges removed
    sum_est_edges[[summary_col_name]] <- rlang::expr(sum(!!rlang::sym(col_name), na.rm = TRUE))
  }
  # Total edges
  sum_est_edges[["n_edges"]] <- rlang::expr(n())
  # Calculate percentage of edges kept: 1 - (total expected edges removed / total edges)
  mut_pct_kept <- list()
  for (col_name in paste0("p_", fracs)) {
    col_name_tot <- paste0(col_name, "_tot")
    col_name_pctkept <- paste0(col_name, "_pctkept")
    mut_pct_kept[[col_name_pctkept]] <- rlang::expr(1 - !!rlang::sym(col_name_tot) / n_edges)
  }
  result_df <- tbl(db$.__enclos_env__$private$con, "edgelist") %>%
    group_by(component) %>%
    mutate(!!!mut_actual_p) %>%
    summarize(!!!sum_est_edges) %>%
    mutate(!!!mut_pct_kept) %>%
    compute(name = "kept_edges", overwrite = TRUE)

  # Compute estimated nodes kept: 1 - (1 - p) ^ degree
  # where 1 - p is the probability of being removed
  sum_nodes <- list()
  for (col_name in paste0("p_", fracs)) {
    col_name_pctkept <- paste0(col_name, "_pctkept")
    col_name_nodes <- paste0(col_name, "_nodes")
    sum_nodes[[col_name_nodes]] <- rlang::expr(sum(1 - ((1 - !!rlang::sym(col_name_pctkept))^degree)))
    sum_nodes[[col_name_pctkept]] <- rlang::expr(mean(!!rlang::sym(col_name_pctkept)))
  }

  # Join the estimated edges kept with the nodes table
  # and estimate the nodes kept
  downsampled_features <- tbl(db$.__enclos_env__$private$con, "nodes") %>%
    left_join(tbl(db$.__enclos_env__$private$con, "kept_edges"), by = "component") %>%
    group_by(component) %>%
    summarize(!!!sum_nodes) %>%
    left_join(node_saturation, by = "component") %>%
    left_join(edge_saturation, by = "component")

  if (verbose && check_global_verbosity()) {
    cli::cli_alert_info("Computing edge and node saturation and average degree for downsampled graphs... ")
  }

  # Compute saturation values and degree
  # Edge saturation: downsampled edges / theoretical max edges
  # Node saturation: downsampled nodes / theoretical max nodes
  # Degree: 2 * downsampled edges / downsampled nodes
  mut_sat <- list()
  for (col_name in paste0("p_", fracs)) {
    col_name_nodes <- paste0(col_name, "_nodes")
    col_name_nodesat <- paste0(col_name, "_nodesat")
    mut_sat[[col_name_nodesat]] <- rlang::expr(!!rlang::sym(col_name_nodes) / theoretical_max_nodes)
    col_name_pct <- paste0(col_name, "_pctkept")
    col_name_edgesat <- paste0(col_name, "_edgesat")
    mut_sat[[col_name_edgesat]] <- rlang::expr((!!rlang::sym(col_name_pct) * as.integer(edges)) / theoretical_max_edges)
    col_name_degree <- paste0(col_name, "_degree")
    mut_sat[[col_name_degree]] <- rlang::expr(2 * (!!rlang::sym(col_name_pct) * as.integer(edges)) / !!rlang::sym(col_name_nodes))
  }

  # Calculate average reads per component to use for
  # converting downsampled proportions to downsampled average reads
  average_reads <- tbl(db$.__enclos_env__$private$con, "edgelist") %>%
    group_by(component) %>%
    summarize(reads = sum(read_count, na.rm = TRUE)) %>%
    pull(reads) %>%
    mean()

  # Apply saturation/degree calculations and format results
  downsampled_saturation <- downsampled_features %>%
    mutate(!!!mut_sat) %>%
    collect() %>%
    select(component, matches("nodesat"), matches("edgesat"), matches("degree")) %>%
    tidyr::pivot_longer(matches("^p_")) %>%
    tidyr::separate(name, into = c("lbl", "p", "type"), sep = "_") %>%
    select(-lbl) %>%
    mutate(p = as.numeric(p)) %>%
    tidyr::pivot_wider(names_from = type, values_from = value) %>%
    mutate(average_reads = round(average_reads * p))

  if (detailed) {
    downsampled_saturation <- list(
      "saturation_downsampled" = downsampled_saturation,
      "saturation_actual" = downsampled_features %>%
        select(component, nodes, node_saturation, theoretical_max_nodes,
               edges, edge_saturation, theoretical_max_edges) %>%
        collect()
    )
  }

  if (verbose && check_global_verbosity()) {
    cli::cli_alert_success("Finished!")
  }

  return(downsampled_saturation)
}
