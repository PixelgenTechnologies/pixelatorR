# Remove CellGraphs

Clears the [`CellGraph`](CellGraph-class.md) objects from the
"cellgraphs" slot.

## Usage

``` r
RemoveCellGraphs(object, ...)

# S3 method for class 'MPXAssay'
RemoveCellGraphs(object, ...)

# S3 method for class 'PNAAssay'
RemoveCellGraphs(object, ...)

# S3 method for class 'PNAAssay5'
RemoveCellGraphs(object, ...)

# S3 method for class 'Seurat'
RemoveCellGraphs(object, assay = NULL, ...)
```

## Arguments

- object:

  An object with cell graphs

- ...:

  Additional arguments (not used)

- assay:

  The name of the target assay

## Value

An object where the elements of the cellgraphs list are set to `NULL`
