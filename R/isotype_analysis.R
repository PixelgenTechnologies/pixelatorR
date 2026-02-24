
#' Residualize a matrix against a model matrix
#'
#' This function takes a matrix M and a model matrix, and returns the residuals
#' of M after regressing out the effects of the model matrix.
#'
#' @param M A matrix of data (cells × features).
#' @param model_mat A matrix of covariates (cells × covariates)
#' @return A matrix of residuals with the same dimensions as M.
#' @noRd
.residualize_matrix <- function(M, model_mat) {
  Q <- qr(model_mat)
  M - model_mat %*% qr.coef(Q, M)
}

#' Partial Least Squares Regression for background factor regression
#'
#' This function performs Partial Least Squares (PLS) regression to model and extract
#' a background score using isotype controls from single-cell data. It can optionally
#' residualize the data against covariates before fitting the PLS model.
#'
#' The PLS model is fit using the isotype control markers as the response variable and
#' the non-isotype markers as predictors. A single factor is extracted from the model,
#' which represents the background signal of all markers in each cell, captured by
#' the markers' covariance with isotype controls. The resulting background score
#' describes the amount of background in each cell, and can be used for downstream
#' analyses or filtering. The model loadings can also be examined to identify which
#' markers are most associated with the background signal. The scores and loadings are
#' oriented such that higher scores correspond to higher background signal, and
#' positive loadings indicate markers that covary positively with the isotype controls.
#'
#' Note that this method assumes that the isotype control markers capture the background
#' signal of interest. If there are other sources of variation that are not captured by
#' the isotype controls, or if the isotype controls do not adequately represent the
#' background signal, then the resulting scores may not fully capture the background
#' variation. If \code{remove_covariates} is \code{TRUE}, the method residualizes the
#' abundance data against the specified covariates before fitting the PLS model. This
#' can help to isolate the background signal by removing variation associated with
#' known covariates (e.g., different samples, cell types, or conditions).
#'
#' @param object Seurat object containing single-cell data.
#' @param isotype_markers A character vector of isotype control marker names.
#' @param model_mat An optional matrix of covariates for residualization (cells × covariates).
#' @param remove_covariates Logical; if TRUE, residualizes data against model_mat before PLS.
#' @param layer String; the data layer in object to use (default is "scale.data").
#'
#' @return A list containing the PLS model, scores, and loadings.
#'
#' @examples
#'
#' library(pixelatorR)
#'
#' seur <-
#' ReadPNA_Seurat(
#'   minimal_pna_pxl_file(),
#'   overwrite = TRUE,
#'   load_proximity_scores = FALSE,
#'   verbose = FALSE
#'   )
#'
#' pls_results <-
#'  isotype_pls(
#'   object = seur,
#'   isotype_markers = c("mIgG1", "mIgG2a", "mIgG2b"),
#'   model_mat = NULL,
#'   remove_covariates = FALSE,
#'   layer = "counts"
#'   )
#'
#' seur$sample <- c("S1", "S1", "S2", "S2", "S2")
#'
#' # Residualizing
#' model_mat <- model.matrix(~ 0 + seur$sample)
#'
#' pls_results_resid <-
#' isotype_pls(
#'  object = seur,
#'  isotype_markers = c("mIgG1", "mIgG2a", "mIgG2b"),
#'  model_mat = model_mat,
#'  remove_covariates = TRUE,
#'  layer = "counts"
#'  )
#'
#' @export
#'
isotype_pls <-
  function(object, isotype_markers, model_mat = NULL, remove_covariates = FALSE, layer = "scale.data") {
    expect_pls()

    assert_class(object, "Seurat")
    assert_vector(isotype_markers, type = "character", n = 1)
    assert_x_in_y(isotype_markers, rownames(object))
    assert_class(model_mat, "matrix", allow_null = TRUE)

    if (!is.null(model_mat) && isTRUE(remove_covariates)) {
      if (nrow(model_mat) != ncol(object)) {
        cli::cli_abort("Number of rows in model_mat must match the number of cells in the Seurat object.")
      }
    }

    assert_single_value(remove_covariates, type = "bool")
    assert_single_value(layer, type = "string")
    assert_x_in_y(layer, Layers(object))

    if (remove_covariates && is.null(model_mat)) {
      cli::cli_abort("model_mat must be provided when remove_covariates is TRUE")
    }

    # transpose to cells × features
    X <-
      object |>
      LayerData(layer = layer) |>
      as.matrix() |>
      t()

    if(nrow(X) == 0 || ncol(X) == 0) {
      cli::cli_abort(glue::glue("No data found in layer '{layer}'. Please check that the layer exists and contains data."))
    }

    # Residualize X
    if (remove_covariates) X <- .residualize_matrix(X, model_mat)

    X_no_isotype <-
      X[, !colnames(X) %in% isotype_markers]
    X_isotype <-
      X[, colnames(X) %in% isotype_markers]

    if (ncol(X_no_isotype) == 0) {
      cli::cli_abort("No non-isotype markers found. At least one non-isotype marker is required for PLS analysis.")
    }

    model <-
      pls::plsr(
        X_isotype ~ X_no_isotype,
        ncomp = 1,
        scale = FALSE,
        validation = "none",
        segments = 10
      )

    score <-
      pls::scores(model)[, 1]

    loadings <-
      pls::loadings(model)[, 1]

    if (cor(score, rowMeans(X_isotype)) < 0) {
      score <- -score
      loadings <- -loadings
    }

    return(
      list(
        model = model,
        scores = score,
        loadings = loadings
      )
    )

  }
