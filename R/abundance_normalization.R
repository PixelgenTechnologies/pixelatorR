# Declarations used in package check
# globalVariables(
#   names = c("modality", "mixture_component"),
#   package = 'pixelatorR',
#   add = TRUE
# )

#' dsb normalization
#'
#' Normalize an MPX count matrix using dsb (\href{https://doi.org/10.1038/s41467-022-29356-8}{Mulè et al, 2022}).
#'
#' @section References:
#' Mulè, M.P., Martins, A.J. & Tsang, J.S. Normalizing and denoising protein expression
#' data from droplet-based single cell profiling. Nat Commun 13, 2099 (2022).
#' https://doi.org/10.1038/s41467-022-29356-8
#'
#' @param object A matrix of MPX counts
#' @param isotype_controls A character vector of isotype controls
#' @param ... Additional arguments. Currently not used.
#'
#' @importFrom stats prcomp
#'
#' @return A matrix of normalized MPX counts
#'
NormalizationMethod.dsb <- function(
    object,
    isotype_controls,
    ...
) {

 if (!requireNamespace(c("mclust", "limma"), quietly = TRUE)) {
    stop(
      "Packages limma and mclust must be installed to use this function.",
      call. = FALSE
    )
  }

  stopifnot(
    "object must be a matrix or Seurat object" =
      inherits(object, c("Matrix", "matrix")),
    "isotype_controls must have at least one element" =
      length(isotype_controls) > 0,
    "isotype_controls must be a character vector" =
      is.character(isotype_controls),
    "All isotype controls must be present in the object" =
      all(isotype_controls %in% rownames(object))
  )

  markers = rownames(object)

  # Get protein negative population means
  object_log = log1p(object)

  protein_model <-
    apply(object_log, 1, function(x) {
      BIC <- mclust::mclustBIC(data = x, G = 2, verbose = FALSE, warn = FALSE)
      mclust::summaryMclustBIC(BIC, x, G = 2)
    })

  mu1 <-
    unlist(lapply(protein_model,
                  function(x) x$parameters$mean[[1]]))

  if (length(mu1) != length(rownames(object_log))) {

    failed_markers <- setdiff(rownames(object_log), names(mu1))

    cli_alert_warning("Empirical background cound not be fit for {length(failed_markers)} proteins: {paste(failed_markers, collapse = ', ')}.")
    cli_alert_info("Values returned will be log transformed without background correction.")

    ad = as.numeric(rep(x = 0, length(failed_markers)))
    names(ad) = failed_markers
    mu1 = c(mu1, ad)
    mu1 = mu1[match(rownames(object_log) , names(mu1) )]
  }

  # Center data to negative mean
  object_norm = apply(object_log, 2, function(x) (x - mu1))

  # Adjust for noise
  cellwise_background_mean <-
    apply(object_norm, 2, function(x) {
      BIC <- mclust::mclustBIC(data = x, G = 2, verbose = FALSE, warn = FALSE)
      return(mclust::summaryMclustBIC(BIC, x, G = 2)$parameters$mean[1])
    })

  noise_matrix <- rbind(object_norm[isotype_controls,],
                       cellwise_background_mean)

  noise_vector <-
    prcomp(t(noise_matrix), scale = TRUE)$x[, 1]

  object_norm <-
    limma::removeBatchEffect(object_norm, covariates = noise_vector)

  return(object_norm)

}

##' CLR normalization
#'
#' Normalize an MPX count matrix using CLR.
#'
#' @param object A matrix of MPX counts
#' @param ... Additional arguments. Currently not used.
#'
#' @return A matrix of normalized MPX counts
#'
NormalizationMethod.clr <- function(
    object,
    ...
) {

  stopifnot(
    "object must be a matrix or Seurat object" =
      inherits(object, c("Matrix", "matrix"))
  )

  sweep(log1p(object), 2, Matrix::colMeans(log1p(object)), "-")

}

#' @rdname NormalizeMPX
#' @method NormalizeMPX Matrix
#'
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

  stopifnot(
    "object must be a matrix or Seurat object" =
      inherits(object, c("Matrix", "matrix"))
  )

  normalization_function <-
    c("dsb" = NormalizationMethod.dsb,
      "clr" = NormalizationMethod.clr)[[method]]

  # Normalize data
  norm_object <-
    normalization_function(object, isotype_controls)

  return(norm_object)
}




#' @rdname NormalizeMPX
#' @method NormalizeMPX MPXAssay
#'
#' @importFrom SeuratObject LayerData LayerData<-
#'
#' @export
#'
NormalizeMPX.MPXAssay <- function(
    object,
    method = c("dsb", "clr"),
    isotype_controls = c("mIgG1", "mIgG2a", "mIgG2b"),
    ...
) {

  newData <-
    NormalizeMPX(object = LayerData(object, "counts"),
                 method = method,
                 isotype_controls = isotype_controls,
                 ...)

  LayerData(object, "data") <-
    as(newData, "CsparseMatrix")

  return(object)

}

#' @rdname NormalizeMPX
#' @method NormalizeMPX CellGraphAssay
#' @docType methods
#' @export
#'
NormalizeMPX.CellGraphAssay <- NormalizeMPX.MPXAssay

#' @rdname NormalizeMPX
#' @method NormalizeMPX CellGraphAssay5
#' @docType methods
#' @export
#'
NormalizeMPX.CellGraphAssay5 <- NormalizeMPX.MPXAssay

#' @rdname NormalizeMPX
#' @method NormalizeMPX Seurat
#'
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
    NormalizeMPX(object = object[[assay]],
                 method = method,
                 isotype_controls = isotype_controls,
                 ...)

  return(object)

}