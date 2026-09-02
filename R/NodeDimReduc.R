#' @include generics.R
NULL

# -------------------------------------------------------
# Class definition
# -------------------------------------------------------

#' The NodeDimReduc class
#'
#' The \code{NodeDimReduc} class stores a dimensionality reduction computed
#' on the nodes of a \code{\link{CellGraph}}. It follows the same design as
#' \code{\link[SeuratObject]{DimReduc}}: a required embeddings matrix, optional
#' feature loadings, standard deviations, a dimension key, and a miscellaneous
#' list for method-specific extras.
#'
#' @slot embeddings A numeric \code{matrix} of node embeddings (nodes x dimensions).
#' Row names are node names and column names use the reduction \code{key}.
#' @slot loadings An optional numeric \code{matrix} of feature loadings
#' (features x dimensions)
#' @slot stdev A numeric vector of standard deviations (or eigenvalues) for
#' each dimension
#' @slot key A character scalar used as the column-name prefix, ending with
#' \code{_} (for example \code{"PC_"})
#' @slot method A character scalar naming the reduction method (for example
#' \code{"pca"} or \code{"umap"})
#' @slot misc A named list of unstructured additional data
#'
#' @name NodeDimReduc-class
#' @rdname NodeDimReduc-class
#' @exportClass NodeDimReduc
#' @concept cellgraph
#'
NodeDimReduc <- setClass(
  Class = "NodeDimReduc",
  slots = list(
    embeddings = "matrix",
    loadings = "matrix",
    stdev = "numeric",
    key = "character",
    method = "character",
    misc = "list"
  ),
  prototype = list(
    embeddings = matrix(nrow = 0, ncol = 0),
    loadings = matrix(nrow = 0, ncol = 0),
    stdev = numeric(),
    key = "DR_",
    method = character(),
    misc = list()
  )
)

# -------------------------------------------------------
# Create methods
# -------------------------------------------------------

#' Create a NodeDimReduc object
#'
#' Constructs a \code{\link{NodeDimReduc}} for storage in the
#' \code{reductions} slot of a \code{\link{CellGraph}}. Analogous to
#' \code{\link[SeuratObject]{CreateDimReducObject}}.
#'
#' @param embeddings A numeric matrix of node embeddings. Row names must be
#' node names. If \code{colnames} are missing they are set to
#' \code{paste0(key, seq_len(ncol(embeddings)))}.
#' @param loadings An optional numeric matrix of feature loadings. Columns are
#' aligned to the embedding dimensions.
#' @param stdev An optional numeric vector of standard deviations, one per
#' dimension
#' @param key Prefix for embedding column names. A trailing underscore is
#' added if missing.
#' @param method Optional name of the reduction method
#' @param misc A list of additional metadata to store with the reduction
#'
#' @return A \code{\link{NodeDimReduc}} object
#'
#' @export
#' @concept cellgraph
#'
CreateNodeDimReducObject <- function(
  embeddings,
  loadings = NULL,
  stdev = numeric(),
  key = "DR_",
  method = character(),
  misc = list()
) {
  assert_class(embeddings, classes = c("matrix", "Matrix"))
  embeddings <- as.matrix(embeddings)
  if (is.null(rownames(embeddings))) {
    cli::cli_abort(c("x" = "{.arg embeddings} must have row names (node names)."))
  }
  if (any(rownames(embeddings) == "") || anyDuplicated(rownames(embeddings))) {
    cli::cli_abort(c("x" = "{.arg embeddings} row names must be unique and non-empty."))
  }

  assert_single_value(key, type = "string")
  if (!grepl("_$", key)) {
    key <- paste0(key, "_")
  }

  if (is.null(colnames(embeddings))) {
    colnames(embeddings) <- paste0(key, seq_len(ncol(embeddings)))
  } else {
    digits <- gsub(paste0("^", key), "", colnames(embeddings))
    if (!all(grepl("^[0-9]+$", digits))) {
      colnames(embeddings) <- paste0(key, seq_len(ncol(embeddings)))
    } else if (!all(startsWith(colnames(embeddings), key))) {
      colnames(embeddings) <- paste0(key, digits)
    }
  }

  if (is.null(loadings) || length(loadings) == 0) {
    loadings <- matrix(nrow = 0, ncol = ncol(embeddings))
    colnames(loadings) <- colnames(embeddings)
  } else {
    assert_class(loadings, classes = c("matrix", "Matrix"))
    loadings <- as.matrix(loadings)
    if (ncol(loadings) != ncol(embeddings)) {
      cli::cli_abort(
        c(
          "x" = "{.arg loadings} must have the same number of columns as {.arg embeddings}.",
          "i" = "Got {ncol(loadings)} and {ncol(embeddings)} columns."
        )
      )
    }
    colnames(loadings) <- colnames(embeddings)
  }

  if (length(stdev) > 0) {
    assert_vector(stdev, type = "numeric", n = 1)
    if (length(stdev) != ncol(embeddings)) {
      cli::cli_abort(
        c(
          "x" = "{.arg stdev} must have one value per embedding dimension.",
          "i" = "Got {length(stdev)} values for {ncol(embeddings)} dimensions."
        )
      )
    }
  }

  new(
    Class = "NodeDimReduc",
    embeddings = embeddings,
    loadings = loadings,
    stdev = as.numeric(stdev),
    key = key,
    method = as.character(method),
    misc = misc %||% list()
  )
}

# -------------------------------------------------------
# Methods
# -------------------------------------------------------

#' NodeDimReduc Methods
#'
#' Methods for \code{\link{NodeDimReduc}} objects
#'
#' @param object A \code{\link{NodeDimReduc}} object
#' @param x A \code{\link{NodeDimReduc}} object
#' @param projected Ignored; included for compatibility with
#' \code{\link[SeuratObject]{Loadings}}
#' @param ... Currently not used
#'
#' @name NodeDimReduc-methods
#' @rdname NodeDimReduc-methods
#' @concept cellgraph
#'
NULL

#' @describeIn NodeDimReduc-methods Show a \code{NodeDimReduc} object
#' @method show NodeDimReduc
#' @docType methods
#'
setMethod(
  f = "show",
  signature = "NodeDimReduc",
  definition = function(object) {
    n_nodes <- nrow(slot(object, "embeddings"))
    n_dims <- ncol(slot(object, "embeddings"))
    method <- slot(object, "method")
    method_lab <- if (length(method) && nzchar(method[1])) method[1] else "dimensional reduction"
    cat(
      "A NodeDimReduc instance of", col_br_blue(method_lab),
      "with", col_br_blue(n_nodes), "nodes and",
      col_br_blue(n_dims), "dimensions\n"
    )
    if (nrow(slot(object, "loadings")) > 0) {
      cat("Feature loadings:", col_br_blue(nrow(slot(object, "loadings"))), "features\n")
    }
  }
)

#' @rdname NodeDimReduc-methods
#' @method Embeddings NodeDimReduc
#' @export
#'
Embeddings.NodeDimReduc <- function(object, ...) {
  slot(object, "embeddings")
}

#' @rdname NodeDimReduc-methods
#' @method Loadings NodeDimReduc
#' @export
#'
Loadings.NodeDimReduc <- function(object, projected = FALSE, ...) {
  slot(object, "loadings")
}

#' @rdname NodeDimReduc-methods
#' @method Stdev NodeDimReduc
#' @export
#'
Stdev.NodeDimReduc <- function(object, ...) {
  slot(object, "stdev")
}

#' @rdname NodeDimReduc-methods
#' @method Key NodeDimReduc
#' @export
#'
Key.NodeDimReduc <- function(object, ...) {
  slot(object, "key")
}

#' @rdname NodeDimReduc-methods
#' @method Cells NodeDimReduc
#' @export
#'
Cells.NodeDimReduc <- function(x, ...) {
  rownames(slot(x, "embeddings"))
}

#' Subset a NodeDimReduc to a set of nodes
#'
#' @noRd
#'
.subset_node_dimreduc <- function(object, nodes) {
  embeddings <- slot(object, "embeddings")
  missing_nodes <- setdiff(nodes, rownames(embeddings))
  if (length(missing_nodes) > 0) {
    cli::cli_abort(
      c(
        "x" = "The reduction is missing embeddings for {length(missing_nodes)} node{?s}.",
        "i" = "Example: {.val {head(missing_nodes, 3)}}"
      )
    )
  }
  slot(object, "embeddings") <- embeddings[nodes, , drop = FALSE]
  object
}

#' Align a NodeDimReduc to a node name order
#'
#' @noRd
#'
.align_node_dimreduc <- function(object, node_names, arg = "reduction") {
  assert_class(object, "NodeDimReduc")
  embeddings <- slot(object, "embeddings")
  slot(object, "embeddings") <- .align_matrix_rows(
    mat = embeddings,
    node_names = node_names,
    arg = arg
  )
  object
}
