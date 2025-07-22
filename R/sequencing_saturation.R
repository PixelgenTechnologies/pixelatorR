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
        res <-
          res %>%
          bind_rows(
            tibble(
              sample_size = sample_sizes[i],
              sample_frac = sample_fracs[i],
              graph_edges = nrow(el),
              graph_proteins = n_distinct(c(el$umi1, el$umi2)),
              graph_reads = length(index)
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
