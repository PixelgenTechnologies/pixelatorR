#' dsb normalization
#'
#' Normalize an MPX count matrix using dsb (\href{https://doi.org/10.1038/s41467-022-29356-8}{Mulè et al, 2022}).
#'
#' @references
#' Mulè, M.P., Martins, A.J. & Tsang, J.S. Normalizing and denoising protein expression
#' data from droplet-based single cell profiling. Nat Commun 13, 2099 (2022).
#'  \url{https://doi.org/10.1038/s41467-022-29356-8}
#'
#' @param counts A matrix of MPX counts
#' @param isotype_controls A character vector of isotype controls
#' @param ... Additional arguments. Currently not used.
#'
#' @examples
#'
#' library(pixelatorR)
#'
#' pxl_file <- minimal_mpx_pxl_file()
#' se <- ReadMPX_Seurat(pxl_file)
#' norm_data_dsb <- NormalizeMPX(
#'   se,
#'   method = "dsb",
#'   isotype_controls = c("mIgG1", "mIgG2a", "mIgG2b")
#' )
#'
#' @return A matrix of normalized MPX counts
#'
#' @noRd
.normalize_method_dsb <- function(
  counts,
  isotype_controls,
  ...
) {
  # Check packages
  expect_mclust()
  expect_limma()

  # Validate input parameters
  assert_class(counts, c("Matrix", "matrix"))
  assert_vector(isotype_controls, type = "character", n = 1)
  assert_x_in_y(isotype_controls, rownames(counts))

  # Get protein negative population means
  counts_log <- log1p(counts)

  protein_model <-
    apply(counts_log, 1, function(x) {
      BIC <- mclust::mclustBIC(data = x, G = 2, verbose = FALSE, warn = FALSE)
      mclust::summaryMclustBIC(BIC, x, G = 2)
    })

  mu1 <-
    unlist(lapply(
      protein_model,
      function(x) x$parameters$mean[[1]]
    ))

  if (length(mu1) != length(rownames(counts_log))) {
    failed_markers <- setdiff(rownames(counts_log), names(mu1))

    cli_alert_warning(
      glue(
        "Empirical background cound not be fit for ",
        "{length(failed_markers)} proteins: ",
        "{paste(failed_markers, collapse = ', ')}."
      )
    )
    cli_alert_info(
      "Values returned will be log transformed without background correction."
    )

    ad <- as.numeric(rep(x = 0, length(failed_markers)))
    names(ad) <- failed_markers
    mu1 <- c(mu1, ad)
    mu1 <- mu1[match(rownames(counts_log), names(mu1))]
  }

  # Center data to negative mean
  counts_norm <- apply(counts_log, 2, function(x) (x - mu1))

  # Adjust for noise
  cellwise_background_mean <-
    apply(counts_norm, 2, function(x) {
      BIC <- mclust::mclustBIC(data = x, G = 2, verbose = FALSE, warn = FALSE)
      return(mclust::summaryMclustBIC(BIC, x, G = 2)$parameters$mean[1])
    })

  noise_matrix <- rbind(
    counts_norm[isotype_controls, ],
    cellwise_background_mean
  )

  noise_vector <-
    prcomp(t(noise_matrix), scale = TRUE)$x[, 1]

  counts_norm <-
    limma::removeBatchEffect(counts_norm, covariates = noise_vector)

  return(counts_norm)
}

##' CLR normalization
#'
#' Normalize an MPX count matrix using CLR.
#'
#' @param counts A matrix of MPX counts
#' @param ... Additional arguments. Currently not used.
#'
#' @examples
#'
#' library(pixelatorR)
#'
#' pxl_file <- minimal_mpx_pxl_file()
#' se <- ReadMPX_Seurat(pxl_file)
#' norm_data_clr <- NormalizeMPX(se, method = "clr")
#'
#' @return A matrix of normalized MPX counts
#'
#' @noRd
.normalize_method_clr <- function(
  counts,
  ...
) {
  assert_class(counts, c("Matrix", "matrix"))

  sweep(log1p(counts), 2, Matrix::colMeans(log1p(counts)), "-")
}

#' @rdname NormalizeMPX
#' @method NormalizeMPX Matrix
#' @docType methods
#' @export
#'
NormalizeMPX.Matrix <- function(
  object,
  method = c("dsb", "clr"),
  isotype_controls = c("mIgG1", "mIgG2a", "mIgG2b"),
  ...
) {
  # Validate inputs
  method <- match.arg(method, choices = c("dsb", "clr"))

  assert_class(object, c("Matrix", "matrix"))

  normalization_function <- switch(
    EXPR = method,
    "dsb" = .normalize_method_dsb,
    "clr" = .normalize_method_clr
  )

  # Normalize data
  norm_object <-
    normalization_function(counts = object, isotype_controls)

  return(norm_object)
}


#' @rdname NormalizeMPX
#' @method NormalizeMPX MPXAssay
#' @docType methods
#' @export
#'
NormalizeMPX.MPXAssay <- function(
  object,
  method = c("dsb", "clr"),
  isotype_controls = c("mIgG1", "mIgG2a", "mIgG2b"),
  ...
) {
  # If object has been merged before, join layers
  if (inherits(object, c("CellGraphAssay5", "Assay5"))) {
    object <- JoinLayers(object)
  }

  newData <-
    NormalizeMPX(
      object = LayerData(object, "counts"),
      method = method,
      isotype_controls = isotype_controls,
      ...
    )

  LayerData(object, "data") <-
    as(newData, "CsparseMatrix")

  return(object)
}

#' @rdname NormalizeMPX
#' @method NormalizeMPX Assay
#' @docType methods
#' @export
#'
NormalizeMPX.Assay <- NormalizeMPX.MPXAssay

#' @rdname NormalizeMPX
#' @method NormalizeMPX CellGraphAssay
#' @docType methods
#' @export
#'
NormalizeMPX.CellGraphAssay <- NormalizeMPX.MPXAssay

#' @rdname NormalizeMPX
#' @method NormalizeMPX Assay5
#' @docType methods
#' @export
#'
NormalizeMPX.Assay5 <- NormalizeMPX.MPXAssay

#' @rdname NormalizeMPX
#' @method NormalizeMPX CellGraphAssay5
#' @docType methods
#' @export
#'
NormalizeMPX.CellGraphAssay5 <- NormalizeMPX.MPXAssay

#' @rdname NormalizeMPX
#' @method NormalizeMPX PNAAssay
#' @docType methods
#' @export
#'
NormalizeMPX.PNAAssay <- function(
  object,
  method = c("dsb", "clr"),
  isotype_controls = c("mIgG1", "mIgG2a", "mIgG2b"),
  ...
) {
  # If object has been merged before, join layers
  if (inherits(object, c("PNAAssay", "PNAAssay5"))) {
    object <- JoinLayers(object)
  }

  newData <-
    NormalizeMPX(
      object = LayerData(object, "counts"),
      method = method,
      isotype_controls = isotype_controls,
      ...
    )

  LayerData(object, "data") <-
    as(newData, "CsparseMatrix")

  return(object)
}

#' @rdname NormalizeMPX
#' @method NormalizeMPX PNAAssay
#' @docType methods
#' @export
#'
NormalizeMPX.PNAAssay <- NormalizeMPX.PNAAssay

#' @rdname NormalizeMPX
#' @method NormalizeMPX Seurat
#' @docType methods
#' @export
#'
NormalizeMPX.Seurat <- function(
  object,
  method = c("dsb", "clr"),
  isotype_controls = c("mIgG1", "mIgG2a", "mIgG2b"),
  assay = NULL,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  object[[assay]] <-
    NormalizeMPX(
      object = object[[assay]],
      method = method,
      isotype_controls = isotype_controls,
      ...
    )

  return(object)
}
