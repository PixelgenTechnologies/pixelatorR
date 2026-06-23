# Normalize MPX data

**\[superseded\]**

This function has been replaced by [`Normalize`](Normalize.md).

## Usage

``` r
NormalizeMPX(
  object,
  method = c("dsb", "clr"),
  isotype_controls = c("mIgG1", "mIgG2a", "mIgG2b"),
  assay = NULL,
  ...
)

# S3 method for class 'Matrix'
NormalizeMPX(
  object,
  method = c("dsb", "clr"),
  isotype_controls = c("mIgG1", "mIgG2a", "mIgG2b"),
  ...
)

# S3 method for class 'MPXAssay'
NormalizeMPX(
  object,
  method = c("dsb", "clr"),
  isotype_controls = c("mIgG1", "mIgG2a", "mIgG2b"),
  ...
)

# S3 method for class 'Assay'
NormalizeMPX(
  object,
  method = c("dsb", "clr"),
  isotype_controls = c("mIgG1", "mIgG2a", "mIgG2b"),
  ...
)

# S3 method for class 'CellGraphAssay'
NormalizeMPX(
  object,
  method = c("dsb", "clr"),
  isotype_controls = c("mIgG1", "mIgG2a", "mIgG2b"),
  ...
)

# S3 method for class 'Assay5'
NormalizeMPX(
  object,
  method = c("dsb", "clr"),
  isotype_controls = c("mIgG1", "mIgG2a", "mIgG2b"),
  ...
)

# S3 method for class 'CellGraphAssay5'
NormalizeMPX(
  object,
  method = c("dsb", "clr"),
  isotype_controls = c("mIgG1", "mIgG2a", "mIgG2b"),
  ...
)

# S3 method for class 'PNAAssay'
NormalizeMPX(
  object,
  method = c("dsb", "clr"),
  isotype_controls = c("mIgG1", "mIgG2a", "mIgG2b"),
  ...
)

# S3 method for class 'PNAAssay5'
NormalizeMPX(
  object,
  method = c("dsb", "clr"),
  isotype_controls = c("mIgG1", "mIgG2a", "mIgG2b"),
  ...
)

# S3 method for class 'Seurat'
NormalizeMPX(
  object,
  method = c("dsb", "clr"),
  isotype_controls = c("mIgG1", "mIgG2a", "mIgG2b"),
  assay = NULL,
  ...
)
```

## Arguments

- object:

  An object.

- method:

  The normalization method to use. Can be "dsb" or "clr".

- isotype_controls:

  A vector of isotype controls to use for normalization.

- assay:

  Name of assay to use; defaults to the default assay.

- ...:

  Additional arguments. Currently not used.

## Value

An object with normalized MPX data.
