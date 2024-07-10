#' Weighted PMDS
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' This function is a wrapper around \code{layout_with_pmds} from
#' the \code{graphlayouts} R package, and can be used to compute a pmds
#' layout of a graph with edge weights.
#'
#' @section weights:
#'
#' - _prob_dist_ : see \code{\link{cos_distance_weights}}
#' - _cos_dist_ : see \code{\link{prob_distance_weights}}
#'
#' @param g An \code{igraph} or a \code{tbl_graph} object
#' @param pivots Number of pivots
#' @param dim Desired number of dimensions. Can be 2 or 3
#' @param method Edge weighting method to use for computing the layout.
#' Can be either "prob_dist" or "cos_dist".
#' @param pow  Power to raise the scores to. Default is 3.
#' @param seed Set seed for pivot sampling
#'
#' @return A matrix of 2D or 3D coordinates
#'
#' @examples
#' library(dplyr)
#'
#' pxl_file <- system.file("extdata/five_cells",
#'                         "five_cells.pxl",
#'                         package = "pixelatorR")
#' seur_obj <- ReadMPX_Seurat(pxl_file) %>%
#'   LoadCellGraphs(cells = colnames(.)[1])
#'
#' # compute weighted pMDS layout
#' g <- CellGraphs(seur_obj)[[1]] %>%
#'   CellGraphData("cellgraph")
#' layout <- layout_with_weighted_pmds(g, dim = 3) %>%
#'   as_tibble(.name_repair = ~c("x", "y", "z"))
#' plotly::plot_ly(
#'   layout,
#'   x = ~x,
#'   y = ~y,
#'   z = ~z,
#'   size = 1,
#'   type = "scatter3d",
#'   mode = "markers"
#' )
#'
#' # or using ComputeLayout and Plot3DGraph
#' seur_obj <- seur_obj %>%
#'   # Compute weighted pMDS layout
#'   ComputeLayout(layout_method = "wpmds", dim = 3)
#'
#' # Create 3D plot
#' Plot3DGraph(seur_obj,
#'             layout_method = "wpmds",
#'             cell_id = colnames(seur_obj)[1],
#'             marker = "CD3E")
#'
#'
#' @export
#'
layout_with_weighted_pmds <- function (
  g,
  dim = 2,
  method = c("prob_dist", "cos_dist"),
  pivots = 200,
  pow  = 3,
  seed = 123
) {

  expect_graphlayouts()

  # Validate input
  stopifnot(
    "'g' must be an 'igraph' object" =
      inherits(g, what = "igraph"),
    "'dim' must be wither 2 or 3" =
      dim %in% c(2, 3),
    "'seed' must be an integer" =
      inherits(seed, what = "numeric"),
    "'pow ' must be a positive numeric value of length 1" =
      inherits(pow , what = "numeric") && (length(pow ) == 1) && (pow  > 0)
  )

  method <- match.arg(method, choices = c("prob_dist", "cos_dist"))

  if (method == "prob_dist") {
    g <- prob_distance_weights(g, k = 5)
  }
  if (method == "cos_dist") {
    g <- cos_distance_weights(g, pivots)
  }

  scores <- g %E>% pull(scores)

  if (pow  != 1) {
    scores <- scores^pow
  }

  set.seed(seed)
  xyz <- graphlayouts::layout_with_pmds(g, pivots = pivots, dim = dim, weights = scores)

  return(xyz)
}

#' Calculate cosine distances
#'
#' @param A,B Matrices with identical dimensions
#'
#' @noRd
cos_dist2 <- function(A, B){
  Matrix::rowSums(A * B) / sqrt(Matrix::rowSums(A * A) * Matrix::rowSums(B * B))
}

#' Calculate edge weights for pMDS
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' @name edge_weights_for_pmds
#' @rdname edge-weights-pmds
#'
#' @param g A \code{tbl_graph} object
#' @param pivots Number of pivots
#'
#' @section Cosine weights:
#' The main idea is that nodes in close spatial proximity should have similar
#' distance profiles to all other nodes in the graph.
#' First, we calculate the pairwise distances between all nodes and a set
#' of randomly selected pivot nodes. After a double centering transformation,
#' we can compute the cosine distance \eqn{d(u, v)} between the distance profiles
#' for nodes \eqn{u} and \eqn{v} connected by an edge in the graph.
#'
#' The cosine distance is a measure of similarity between
#' two vectors and ranges from -1 to 1, where 1 indicates maximal similarity and
#' -1 indicates maximal dissimilarity. Since we want to assign higher weights to
#' nodes that are dissimilar, we convert the distances as follows:
#'
#' \deqn{w_{u,v} = 2 - d(u, v)}
#'
#' @export
#'
cos_distance_weights <- function (
  g,
  pivots
) {

  stopifnot(
    "'g' must be a 'tbl_graph' object" =
      inherits(g, what = "tbl_graph"),
    "'pivots' must be a positive numeric value" =
      inherits(pivots, what = "numeric") &&
      (pivots > 0) & (length(pivots) == 1)
  )

  pivs <- sample(1:igraph::vcount(g), pivots)

  D <- t(igraph::distances(g, v = pivs, weights = NA))
  D <- D^2
  cmean <- Matrix::colMeans(D)
  rmean <- Matrix::rowMeans(D)
  Dmat <- D - outer(rmean, cmean, function(x, y) x + y) +
    mean(D)

  # Calculate cosine distances
  cdist2 <- cos_dist2(Dmat[g %E>% pull(from), ], Dmat[g %E>% pull(to), ])

  # fetch scores and place in graph edge table
  scores <- 2 - cdist2
  g <- g %E>%
    select(-any_of(c("scores"))) %>%
    mutate(scores = scores)

  return(g)
}

#' @rdname edge-weights-pmds
#'
#' @param k Number of steps
#'
#' @section Transition probability weights:
#' Weights are calculated for edges in the graph from the bidirectional
#' transition probability of a random walker. For an undirected edge \eqn{e_{i,j}},
#' we define the bidirectional transition probability \eqn{P(e_{u,v}} as the
#' probability of transitioning from \eqn{u} to \eqn{v} and back to \eqn{u}
#' in \eqn{k} steps. Note that it doesn't matter what node the random walker
#' starts from since it has to go in both directions. In this computation,
#' we also allow self-loops in the graph to slow the transition of the
#' random walker.
#'
#' Then, we assume that there is a relationship between the bidirectional
#' transition probability and an exponential function which decays with distance:
#'
#' \deqn{P(e_{u,v}) \propto e^{-d(u, v)}}
#'
#' by solving for \eqn{d(u, v)}, we can define an edge weight that is proportional
#' to distance:
#'
#' \deqn{w_{u,v} = -\log(P(e_{u,v}))}
#'
#' @export
#'
prob_distance_weights <- function (
  g,
  k = 5
) {

  stopifnot(
    "'g' must be a 'tbl_graph' object" =
      inherits(g, what = "tbl_graph")
  )

  A <- as_adjacency_matrix(g)
  diag(A) <- 1
  P <- A / Matrix::rowSums(A)

  P_steps <- Reduce("%*%", rep(list(P), k))
  P_steps <- P_steps * A
  P_steps_bidirectional <- P_steps * Matrix::t(P_steps)

  # Extract edge scores and place in graph edge table
  r_ids <- P_steps_bidirectional@i + 1
  c_ids <- findInterval(seq(P_steps_bidirectional@x) - 1, P_steps_bidirectional@p[-1]) + 1
  map_vals <- tibble(
    from = r_ids,
    to = c_ids,
    bi_prob = P_steps_bidirectional@x
  )
  g <- g %E>%
    select(-any_of(c("bi_prob", "scores"))) %>%
    left_join(y = map_vals, by = c("from", "to"))
  g <- g %E>% mutate(scores = -log(bi_prob))

  return(g)
}

