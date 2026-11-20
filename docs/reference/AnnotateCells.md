# Automatic annotation of cell types

This function transfers cell type annotations from a reference Seurat
dataset to a query Seurat dataset using non-negative least squares
(`"nnls"`, default) or Seurat anchor transfer (`"Seurat"`). For each
column specified in `reference_groups`, a corresponding column of
predicted cell type annotations and a `<group>_score` column of per-cell
prediction scores are added to the metadata of the query dataset. In
addition, one can optionally summarize these annotations for a specific
clustering solution, which adds a `<group>_summary` column containing
the most common annotation per cluster.

## Usage

``` r
AnnotateCells(
  object,
  reference = reference,
  summarize_by_column = NULL,
  reference_assay = "PNA",
  query_assay = "PNA",
  reference_groups = c("celltype.l1", "celltype.l2"),
  normalization_method = c("CLR", "LogNormalize", "SCT"),
  reduction = "cca",
  method = c("nnls", "Seurat", "nmf"),
  min_prediction_score = 0,
  predict_doublets = FALSE,
  min_doublet_score = 0.4,
  skip_normalization = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  Seurat object to which you want to add celltype annotation

- reference:

  Seurat object that contains annotated reference data, see
  [`read_pbmc_reference`](read_pbmc_reference.md) and examples for
  obtaining a PBMC reference dataset

- summarize_by_column:

  If NULL, return only the cell type annotation per cell. If the name of
  a column in the object metadata is provided, the annotations will be
  summarized to the factors of this column, i.e. the most common
  annotation per factor level will be retained. Useful for annotating
  the output of the Seurat clustering, for example, Default: NULL

- reference_assay:

  Assay in the reference dataset used for annotation, Default: 'PNA'

- query_assay:

  Assay in the query dataset used for annotation, Default: 'PNA'

- reference_groups:

  A character vector of reference metadata columns to use for
  annotation, Default: `c("celltype.l1", "celltype.l2")`

- normalization_method:

  Normalization method used prior to annotation (with margin = 2),
  either "CLR", "LogNormalize", or "SCT", Default: 'CLR'

- reduction:

  Reduction method used during `FindTransferAnchors` when
  `method = "Seurat"`, Default: 'cca'

- method:

  Method to use for annotation, either "nnls", "Seurat", or "nmf" (alias
  for "nnls"), Default: 'nnls'

- min_prediction_score:

  Minimum prediction score for the annotation to be considered valid.
  Labels with a prediction score below this value will be labelled as
  "Unknown", Default: 0, meaning that no threshold is used.

- predict_doublets:

  If TRUE, compute a heterotypic doublet score and candidate cell type
  pair using the first group in `reference_groups` when
  `method = "nnls"`, adding `doublet_score` and `doublet_pair` columns
  to the object metadata, Default: FALSE

- min_doublet_score:

  Minimum doublet score required to report a candidate cell type pair in
  `doublet_pair` when `predict_doublets = TRUE`. Cells with a doublet
  score below this threshold will have `doublet_pair` set to `NA`,
  Default: 0.4

- skip_normalization:

  If TRUE, the datasets will not be normalized prior to annotation. The
  default layer will be used, so use this setting at your own risk,
  Default: FALSE

- verbose:

  If TRUE, print messages about the progress of the annotation, Default:
  TRUE

- ...:

  Additional parameters passed to `FindTransferAnchors` when
  `method = "Seurat"`

## Value

The initial Seurat object with additional annotation columns added to
its metadata.

## Details

The "nnls" (or "nmf") method projects query cells onto feature-weighted
reference cell type profiles using non-negative least squares
([`RcppML::project`](https://rdrr.io/pkg/RcppML/man/project.html)),
followed by local fine-tuning across the top candidate cell types. The
"Seurat" method is a wrapper for the `FindTransferAnchors` and
`TransferData` functions from Seurat. Both methods support an optional
summary per cluster.

Merged SeuratObject v5 assays (`Assay5`, including `PNAAssay5`) must not
retain split merge layers (`counts.1`, `data.1`, etc.). Run
`JoinLayers()` on the `Seurat` object (or on the assay) after
[`merge()`](https://rdrr.io/r/base/merge.html) before calling
`AnnotateCells()`.

## Examples

``` r
if (FALSE) { # \dontrun{
reference <- read_pbmc_reference()

seur <- ReadPNA_Seurat(minimal_pna_pxl_file())

seur <- AnnotateCells(
  seur,
  reference = reference
)
} # }
```
