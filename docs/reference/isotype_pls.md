# Partial Least Squares Regression for background factor regression

This function performs Partial Least Squares (PLS) regression to model
and extract a background score using isotype controls from single-cell
data. It can optionally residualize the data against covariates before
fitting the PLS model.

## Usage

``` r
isotype_pls(
  object,
  isotype_markers,
  model_mat = NULL,
  remove_covariates = FALSE,
  layer = "scale.data"
)
```

## Arguments

- object:

  Seurat object containing single-cell data.

- isotype_markers:

  A character vector of isotype control marker names.

- model_mat:

  An optional matrix of covariates for residualization (cells ×
  covariates).

- remove_covariates:

  Logical; if TRUE, residualizes data against model_mat before PLS.

- layer:

  String; the data layer in object to use (default is "scale.data").

## Value

A list containing the PLS model, scores, and loadings.

## Details

The PLS model is fit using the isotype control markers as the response
variable and the non-isotype markers as predictors. A single factor is
extracted from the model, which represents the background signal of all
markers in each cell, captured by the markers' covariance with isotype
controls. The resulting background score describes the amount of
background in each cell, and can be used for downstream analyses or
filtering. The model loadings can also be examined to identify which
markers are most associated with the background signal. The scores and
loadings are oriented such that higher scores correspond to higher
background signal, and positive loadings indicate markers that covary
positively with the isotype controls.

Note that this method assumes that the isotype control markers capture
the background signal of interest. If there are other sources of
variation that are not captured by the isotype controls, or if the
isotype controls do not adequately represent the background signal, then
the resulting scores may not fully capture the background variation. If
`remove_covariates` is `TRUE`, the method residualizes the abundance
data against the specified covariates before fitting the PLS model. This
can help to isolate the background signal by removing variation
associated with known covariates (e.g., different samples, cell types,
or conditions).

## Examples

``` r
library(pixelatorR)

seur <-
  ReadPNA_Seurat(
    minimal_pna_pxl_file(),
    overwrite = TRUE,
    load_proximity_scores = FALSE,
    verbose = FALSE
  )

pls_results <-
  isotype_pls(
    object = seur,
    isotype_markers = c("mIgG1", "mIgG2a", "mIgG2b"),
    model_mat = NULL,
    remove_covariates = FALSE,
    layer = "counts"
  )

seur$sample <- c("S1", "S1", "S2", "S2", "S2")

# Residualizing
model_mat <- model.matrix(~ 0 + seur$sample)

pls_results_resid <-
  isotype_pls(
    object = seur,
    isotype_markers = c("mIgG1", "mIgG2a", "mIgG2b"),
    model_mat = model_mat,
    remove_covariates = TRUE,
    layer = "counts"
  )
```
