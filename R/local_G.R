#' Calculate Local G
#'
#' Local G is a local spatial statistic that measures the degree of clustering protein marker counts.
#'
#' @section Details:
#' Local G can for instance be used to detect hot spots for protein markers in a graph, where nodes
#' that are close to each other have similar marker count values. The metric is a Z-score that
#' measures the deviation of the observed local marker expression from the expected marker
#' expression under the null hypothesis of no spatial association. The sign of the score
#' indicates whether the observed marker counts are higher or lower than expected
#' and can therefore be used to identify hot and/or cold spots.
#'
#' The observed local marker expression for a node is computed by aggregating the weighted marker
#' expression within its local neighborhood. The local G metric is largely influenced by the choice
#' of edge weights (see section on weights below) and the size of the local neighborhood \code{k}.
#' The method implemented here uses incoming transition probabilities for a k-step random walk as
#' edge weights. By increasing \code{k}, the local neighborhood of each node is expanded, increasing
#' the "smoothing effect" of the metric, which can be useful to increase the scale of spatial association.
#'
#' Local G can also be useful for more interpretative visualization of polarized marker expression,
#' as it enhances spatial trends across neighborhoods in the graph, even if the marker counts in
#' individual nodes are sparse.
#'
#' @section Definition of \eqn{G_{i}} ("gi"):
#' Definition of local G (Ord and Getis 1995, equation 6):
#' \deqn{
#' Z(G_i)=\dfrac{[\sum_{j=1}^{n}w_{i,j}x_j]-[(\sum_{j=1}^{n}w_{i,j})\bar x^*]}
#' {s^*\{[(n-1)\sum_{j=1}^{n}w_{i,j}^2-(\sum_{j=1}^{n}w_{i,j})^2]/(n-2)\}^{1/2}},\forall j\neq i
#' }
#' \cr
#' where \deqn{s_i = \sqrt{((\sum_{j=1}^{n}x_j^2)/(n-1))-[\bar x_i]^2},\forall j\neq i}
#' and \deqn{\bar x_i=(\sum_{j=1}^{n}x_j)/(n-1),\forall j\neq i}
#'
#' @section Definition of \eqn{G_{i}^{*}} ("gstari"):
#' In the equation for \eqn{G_{i}}, the condition that \eqn{i\neq j} is central. \eqn{G_i^*}
#' relaxes this constraint, by including \eqn{i} as a neighbor of itself. This statistic
#' is expressed as (Ord and Getis 1995, equation 7):
#' \deqn{
#' Z(G_i^*)=\dfrac{[\sum_{j=1}^{n}w_{i,j}x_j]-[\sum_{j=1}^{n}w_{i,j}\bar x_i]}
#' {s_i\{[n\sum_{j=1}^{n}w_{i,j}^2-(\sum_{j=1}^{n}w_{i,j})^2]/(n-1)\}^{1/2}}
#' }
#' \cr
#' where \deqn{s^* = \sqrt{((\sum_{j=1}^{n}x_j^2)/n)-\bar x_i^{*2}}}
#' and \deqn{\bar x_i^*=(\sum_{j=1}^{n}x_j)/n}
#' The left numerator corresponds to \eqn{G_i}, the right to \eqn{E(G_i)},
#' and the denominator to \eqn{Var(G_i)}.
#'
#' @section Weights:
#' Weights are calculated for node pairs if \code{use_weights = TRUE}. For a center node \eqn{i},
#' we calculate weights for each neighbor \eqn{j \in N(i)} as the the transition probability
#' for a k-step walk from \eqn{j} to \eqn{i}. This helps normalizing the contribution
#' by each neighbor in a manner that reduces the influence of high degree nodes. It also
#' ensures that nodes that are further away from the central node in each neighborhood
#' are given lower weights. The strategy rests on the assumption that node degree correlates
#' strongly with the node marker count.
#'
#' In practice, the weights are calculated by normalizing the adjacency matrix of the graph
#' such that each row sums to 1 (also known as the stochastic matrix). By multiplying the
#' transpose of the stochastic matrix with the count matrix, we obtain a "lag matrix" with the weighted
#' marker expression for each node neighborhood which is leveraged in the computation of local G.
#'
#' When using \code{k > 1}, the \code{k}'th power of the stochastic matrix is used instead. In this
#' case, the transition probabilities are calculated for a random walk of length \code{2} or higher.
#' Since the random walker can reach a larger neighborhood when \code{k > 1}, we will
#' also include information about the marker expression for nodes that are more than 1 step
#' away from the center node. Moreover, when \code{k > 1}, the weights for \eqn{i \sim i}
#' (self-loops) can become positive. These weights are not allowed in the calculation of
#' \eqn{G_i} and are therefore set to 0, and the remaining transition probabilities are normalized
#' to 1 within each node neighborhood.
#'
#' The weighting scheme is ignored if \code{W} is provided.
#'
#' @param g A \code{tbl_graph} object representing an MPX/PNA component graph
#' @param counts A \code{dgCMatrix} (sparse matrix) with node marker counts. Raw counts
#' are recommended.
#' @param W A \code{dgCMatrix} (sparse matrix) with edge weights. If not provided,
#' the function will use the adjacency matrix of the graph to compute weights.
#' @param use_weights Logical indicating whether to use edge weights. See section
#' \strong{weights} for details.
#' @param normalize_counts Logical indicating whether to normalize the counts. If
#' \code{TRUE}, the counts within each node are converted to proportions.
#' @param type One of "gstari" or "gi". See details below for more information.
#' @param return_p_vals Compute and return p-values for the local G scores. The
#' p values are also adjusted for multiple testing per marker.
#' @param p_adjust_method Method to use for multiple testing adjustment. See \code{\link{p.adjust}}
#' for details.
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of "two.sided" (default), "greater" or "less".
#' @param ... Additional arguments passed from other methods.
#'
#' @inheritParams compute_transition_probabilities
#'
#' @references
#' - Ord, J. K., Getis. A. Local Spatial Autocorrelation Statistics:
#' Distributional Issues and an Application
#' \url{https://doi.org/10.1111/j.1538-4632.1995.tb00912.x}
#'
#' @examples
#' library(dplyr)
#' # Load example data as a Seurat object
#' pxl_file <- minimal_mpx_pxl_file()
#' seur_obj <- ReadMPX_Seurat(pxl_file) %>%
#'   LoadCellGraphs(cells = colnames(.)[1]) %>%
#'   ComputeLayout(layout_method = "pmds", dim = 3)
#' cg <- CellGraphs(seur_obj)[[1]]
#' g <- CellGraphData(cg, "cellgraph")
#' counts <- CellGraphData(cg, "counts")
#' xyz <- CellGraphData(cg, "layout")[["pmds_3d"]]
#'
#' # Compute local G scores
#' gi_mat <- local_G(g, counts = counts)
#'
#' # Visualize
#' scores <- gi_mat[, "CD3E"]
#' max_abs_val <- max(abs(scores))
#' node_colors <-
#'   scales::col_quantile(
#'     domain = c(-max_abs_val, max_abs_val),
#'     palette = "RdBu",
#'     probs = seq(0, 0.99, length.out = 11),
#'     reverse = TRUE
#'   )(scores)
#' plotly::plot_ly(
#'   xyz,
#'   x = ~x,
#'   y = ~y,
#'   z = ~z,
#'   type = "scatter3d",
#'   mode = "markers",
#'   marker = list(
#'     size = 5,
#'     color = node_colors
#'   )
#' )
#'
#' @return A matrix with local G scores or a list with the following items:
#' \itemize{
#'    \item{gi_mat: A matrix with local G scores}
#'    \item{gi_p_mat: A matrix with p-values (if \code{return_p_vals = TRUE})}
#'    \item{gi_p_adj_mat: A matrix with adjusted p-values (if \code{return_p_vals = TRUE})}
#' }
#'
#' @export
#'
local_G <- function(
  g,
  counts,
  k = 1,
  W = NULL,
  use_weights = TRUE,
  normalize_counts = FALSE,
  type = c("gi", "gstari"),
  return_p_vals = FALSE,
  p_adjust_method = "BH",
  alternative = c("two.sided", "less", "greater"),
  ...
) {
  # Check input parameters
  assert_class(g, c("tbl_graph", "igraph"))
  assert_class(counts, c("Matrix", "matrix"))
  assert_single_value(k, type = "integer")
  if (!k > 0) {
    cli::cli_abort(
      c(
        "i" = "{.var k} must be a positive {.cls integer}",
        "x" = "k = {k}"
      )
    )
  }
  assert_singles_match(nrow(counts), length(g))
  assert_single_value(use_weights, type = "bool")
  assert_single_value(normalize_counts, type = "bool")
  assert_single_value(return_p_vals, type = "bool")

  type <- match.arg(type, choices = c("gi", "gstari"))
  alternative <- match.arg(alternative, choices = c("two.sided", "less", "greater"))

  # Get adjacency matrix
  A <- as_adjacency_matrix(g)

  # Normalize marker counts
  if (normalize_counts) {
    counts <- sweep(counts, 1, Matrix::rowSums(counts), "/")
  }

  # Validate W
  if (!is.null(W)) {
    assert_class(W, "dgCMatrix")
    assert_singles_match(nrow(W), length(g))
    assert_singles_match(ncol(W), length(g))
    if (type == "gstari" && any(diag(A) == 0)) {
      cli::cli_abort(
        c("x" = "The 'gstari' type requires the diagonal of the {.var A} matrix to be positive.")
      )
    }
  } else {
    if (use_weights) {
      # If gstari is used, allow self-loops
      if (type == "gstari") {
        diag(A) <- 1
      }
      # Compute transition probabilities
      W <- compute_transition_probabilities(A, k = k, remove_self_loops = (type == "gi")) %>% Matrix::t()
    } else {
      # If no normalization is desired, use the adjacency matrix where
      # each weight is 1
      W <- A
      if (k > 1) {
        # Exapand local neighborhood using matrix powers
        W <- Reduce(`%*%`, rep(list(W), k))
        W@x <- rep(1, length(W@x))
        if (type == "gstari") {
          diag(W) <- 1
        } else {
          diag(W) <- 0
        }
      }
    }
  }

  # get number of nodes
  n_nodes <- length(g)

  # Calculate lag matrix
  lag_mat <- W %*% counts

  if (type == "gstari") {
    # Compute xibar for each node and marker
    xibar_mat <- apply(counts, 2, function(val) {
      rep(mean(val), n_nodes)
    })
    # Compute si for each node and marker
    s_mat <- do.call(cbind, lapply(seq_len(ncol(counts)), function(i) {
      rep((sum(counts[, i]^2) / n_nodes) - (sum(counts[, i]) / n_nodes)^2, n_nodes)
    }))
  } else if (type == "gi") {
    # Compute xibar for each node and marker, excluding i
    xibar_mat <- apply(counts, 2, function(val) {
      (rep(sum(val), n_nodes) - val) / (n_nodes - 1)
    })
    # Compute si for each node and marker, excluding i
    s_mat <- do.call(cbind, lapply(seq_len(ncol(counts)), function(i) {
      ((rep(sum(counts[, i]^2), n_nodes) - counts[, i]^2) / (n_nodes - 1)) - xibar_mat[, i]^2
    }))
  }

  # Calculate weights for each node. If k = 1, the weights
  # should be equal to the node degree
  weights_i <- rowSums(W)
  square_weights_i <- rowSums(W^2)

  # Calculate expected value
  E_G <- weights_i * xibar_mat

  # Calculate numerator G - E(G)
  numerator_mat <- (lag_mat - E_G)

  # Calculate denomenator Var(G)
  if (type == "gstari") {
    denomenator_mat <- sqrt(s_mat * ((n_nodes * square_weights_i - weights_i^2) / (n_nodes - 1)))
  } else if (type == "gi") {
    denomenator_mat <- sqrt(s_mat * (((n_nodes - 1) * square_weights_i - weights_i^2) / (n_nodes - 2)))
  }

  # Calculate Z-score
  gi_mat <- as.matrix(numerator_mat / denomenator_mat)

  # Set non finite values to 0
  gi_mat[!is.finite(gi_mat)] <- 0

  # Handle p-value computation
  if (return_p_vals) {
    if (alternative == "two.sided") {
      gi_p_mat <- apply(gi_mat, 2, function(x) {
        2 * pnorm(abs(x), lower.tail = FALSE)
      })
    } else if (alternative == "greater") {
      gi_p_mat <- apply(gi_mat, 2, function(x) {
        pnorm(x, lower.tail = FALSE)
      })
    } else if (alternative == "less") {
      gi_p_mat <- apply(gi_mat, 2, function(x) {
        pnorm(x)
      })
    }
    gi_p_adj_mat <- apply(gi_p_mat, 2, p.adjust, method = p_adjust_method)
    return(list(gi_mat = gi_mat, gi_p_mat = gi_p_mat, gi_p_adj_mat = gi_p_adj_mat))
  }

  return(gi_mat)
}


#' Compute transition probabilities
#'
#' Create a transition probability matrix from an adjacency matrix.
#'
#' @param A An adjacency matrix
#' @param k `r lifecycle::badge("experimental")` An integer value specifying the maximum steps
#' from each node to expand the neighborhood. E.g. with \code{k = 2}, each node neighborhood
#' will include the direct neighbors and nodes two steps away. With \code{k = 1} (default), only
#' the direct neighbors will be used.
#' @param remove_self_loops Whether to remove self-loops from the transition probability matrix.
#'
#' @examples
#' library(igraph)
#' g <- make_lattice(c(2, 3))
#' A <- as_adjacency_matrix(g)
#'
#' # Transition probabilities
#' P_out <- compute_transition_probabilities(A)
#' P_out
#'
#' @return A matrix of transition probabilities. The transition probability \eqn{P^k(u \rightarrow v)}
#' for a k-step walk is found in row \eqn{u} and column \eqn{v} of the transition probability matrix.
#' Row \eqn{v} and column \eqn{u} gives the reversed transition probability \eqn{P^k(v \rightarrow u)}
#' for the k-step walk.
#'
#' @export
#'
compute_transition_probabilities <- function(
  A,
  k = 1,
  remove_self_loops = FALSE
) {
  if (nrow(A) != ncol(A)) {
    cli::cli_abort(
      c(
        "i" = "{.var A} must be a square matrix",
        "x" = "nrow(A) = {nrow(A)}, ncol(A) = {ncol(A)}"
      )
    )
  }
  if (!Matrix::isSymmetric(A)) {
    cli::cli_abort(
      c("x" = "{.var A} must be a symmetric matrix")
    )
  }
  assert_single_value(k, type = "integer")
  if (!k > 0) {
    cli::cli_abort(
      c(
        "i" = "{.var k} must be a positive {.cls integer}",
        "x" = "k = {k}"
      )
    )
  }

  # Compute transition probabilities
  W_out <- A * (1 / Matrix::rowSums(A))
  # Compute matrix power determined by k
  W_out <- Reduce("%*%", rep(list(W_out), k))
  # When using gi, set diagonal to 0 and recompute probabilities
  if (remove_self_loops && (k > 1)) {
    diag(W_out) <- 0
    W_out <- W_out * (1 / Matrix::rowSums(W_out))
  }

  return(W_out)
}
