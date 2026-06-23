# Normalize MPX or PNA data

Normalizes MPX or PNA data using the specified method. The normalization
method can be one of "dsb" or "CLR".

## Usage

``` r
Normalize(
  object,
  method = c("dsb", "clr"),
  isotype_controls = c("mIgG1", "mIgG2a", "mIgG2b"),
  assay = NULL,
  ...
)

# S3 method for class 'Matrix'
Normalize(
  object,
  method = c("dsb", "clr"),
  isotype_controls = c("mIgG1", "mIgG2a", "mIgG2b"),
  ...
)

# S3 method for class 'MPXAssay'
Normalize(
  object,
  method = c("dsb", "clr"),
  isotype_controls = c("mIgG1", "mIgG2a", "mIgG2b"),
  ...
)

# S3 method for class 'Assay'
Normalize(
  object,
  method = c("dsb", "clr"),
  isotype_controls = c("mIgG1", "mIgG2a", "mIgG2b"),
  ...
)

# S3 method for class 'CellGraphAssay'
Normalize(
  object,
  method = c("dsb", "clr"),
  isotype_controls = c("mIgG1", "mIgG2a", "mIgG2b"),
  ...
)

# S3 method for class 'Assay5'
Normalize(
  object,
  method = c("dsb", "clr"),
  isotype_controls = c("mIgG1", "mIgG2a", "mIgG2b"),
  ...
)

# S3 method for class 'CellGraphAssay5'
Normalize(
  object,
  method = c("dsb", "clr"),
  isotype_controls = c("mIgG1", "mIgG2a", "mIgG2b"),
  ...
)

# S3 method for class 'PNAAssay'
Normalize(
  object,
  method = c("dsb", "clr"),
  isotype_controls = c("mIgG1", "mIgG2a", "mIgG2b"),
  ...
)

# S3 method for class 'PNAAssay'
Normalize(
  object,
  method = c("dsb", "clr"),
  isotype_controls = c("mIgG1", "mIgG2a", "mIgG2b"),
  ...
)

# S3 method for class 'Seurat'
Normalize(
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

An object with normalized MPX or PNA data.

## Details

CLR can be used to normalize MPX or PNA data using the centered
log-ratio transformation in which the assumption is that the geometric
mean of the marker abundance is constant across cells (e.g cell lines).
This assumption might not hold for data sets from samples from different
sources or having a variable cell type composition. In addition, CLR
does not take the noise from unspecific binding of antibodies into
account.

For these reasons, the dsb normalization method can be a useful
alternative in mixed-population data sets. dsb normalizes marker counts
based on their abundance in a negative population across all cells and
regresses out a per-cell noise component based on isotype controls and
non-specific marker abundance.

## References

Mulè, M.P., Martins, A.J. & Tsang, J.S. Normalizing and denoising
protein expression data from droplet-based single cell profiling. Nat
Commun 13, 2099 (2022). <https://doi.org/10.1038/s41467-022-29356-8>
