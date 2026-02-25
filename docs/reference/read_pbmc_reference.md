# Read PBMC cell annotation reference

This function reads the cell annotation reference data from the
pixelatorR package's extdata directory. This data contains 720 annotated
PBMC components with their corresponding counts of 158 proteins. The
annotation has been performed at two levels of granularity:
"celltype.l1" (broad categories) and "celltype.l2" (more specific
subtypes). When loading the data through this function, a Seurat object
is created with the counts data stored in the "PNA" assay and the cell
annotations included in the metadata.

## Usage

``` r
read_pbmc_reference()
```

## Value

A Seurat object containing the cell annotation and counts data.
