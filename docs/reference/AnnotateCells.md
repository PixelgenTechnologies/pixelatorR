# Automatic annotation of cell types

This function finds anchors in a reference dataset and transfers the
cell type of the reference to a query Seurat dataset. Two columns of
cell type annotations are added to the metadata of the query dataset,
one high-level and one more fine-grained. In addition, one can
optionally summarize these annotations for a specific clustering
solution. When choosing this, the most common annotation per cluster is
added to the metadata of the object.

## Usage

``` r
AnnotateCells(
  object,
  reference = reference,
  summarize_by_column = NULL,
  reference_assay = "ADT",
  query_assay = "PNA",
  reference_groups = c("celltype.l1", "celltype.l2"),
  normalization_method = c("LogNormalize", "SCT", "CLR"),
  reduction = "cca",
  method = c("Seurat", "nmf"),
  min_prediction_score = 0,
  skip_normalization = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  Seurat object to which you want to add celltype annotation

- reference:

  Seurat object that contains annotated reference data, see examples for
  a programmatic way to obtain a PBMC reference dataset

- summarize_by_column:

  If NULL, return only the cell type annotation per cell. If the name of
  a column in the object metadata is provided, the annotations will be
  summarized to the factors of this column, i.e. the most common
  annotation per factor level will be retained. Useful for annotating
  the output of the Seurat clustering, for example.

- reference_assay:

  Assay in the reference dataset used to find the anchors, Default:
  'ADT'

- query_assay:

  Assay in the query dataset used to find the anchors. It should
  probably be set to a normalized assay that was run on the object (e.g.
  "dsb"), Default: 'PNA'

- reference_groups:

  A character vector of reference groups to use for annotation.

- normalization_method:

  normalization method used during `FindTransferAnchors`, Default:
  'LogNormalize'

- reduction:

  reduction method used during `FindTransferAnchors`, Default: 'cca'

- method:

  Method to use for annotation, either "Seurat" or "nmf".

- min_prediction_score:

  Minimum prediction score for the annotation to be considered valid.
  Labels with a prediction score below this value will be labelled as
  "Unknown". Default: 0, meaning that not threshold is used.

- skip_normalization:

  If TRUE, the datasets will not be normalized prior to annotation. This
  parameter only has an effect if `method` is set to "nmf". The default
  layer will be used, so use this setting at your own risk.

- verbose:

  If TRUE, print messages about the progress of the annotation.

- ...:

  Additional parameters to `FindTransferAnchors`

## Value

The initial Seurat object with additional annotation columns added to
its metadata.

## Details

The "Seurat" method is a wrapper for the `FindTransferAnchors` and
`TransferData` functions from Seurat, followed by an optional summary
per cluster.

## Examples

``` r
if (FALSE) { # \dontrun{
# Download reference file
reference <-
  readRDS(url(
    paste0(
      "https://pixelgen-technologies-datasets.s3.eu-north-1.amazonaws.com/",
      "mpx-analysis/next/R/pbmc_annotation.rds"
    )
  ))

# Does not work on the test data due to small size - ignore for now.
seur <- ReadPNA_Seurat(minimal_pna_pxl_file())

seur <- AnnotateCells(seur,
  reference = reference, query_assay = "PNA"
)
} # }
```
