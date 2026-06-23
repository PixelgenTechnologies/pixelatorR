# CellGraphAssay Methods

Methods for [`CellGraphAssay`](CellGraphAssay-class.md) objects for
generics defined in other packages

## Usage

``` r
# S3 method for class 'CellGraphAssay'
RenameCells(object, new.names = NULL, ...)

# S4 method for class 'CellGraphAssay'
show(object)

# S3 method for class 'CellGraphAssay'
subset(x, features = NULL, cells = NULL, ...)

# S3 method for class 'CellGraphAssay'
merge(
  x = NULL,
  y = NULL,
  merge.data = TRUE,
  add.cell.ids = NULL,
  collapse = TRUE,
  ...
)
```

## Arguments

- object:

  A `CellGraphAssay` or a `CellGraphAssay5` object

- new.names:

  A character vector with new cell IDs. The length of the vector must be
  equal to the number of cells in the object and the names must be
  unique.

- ...:

  Arguments passed to other methods

- x:

  A [`CellGraphAssay`](CellGraphAssay-class.md) object

- features:

  Feature names

- cells:

  Cell names

- y:

  A [`CellGraphAssay`](CellGraphAssay-class.md) object or a list of
  [`CellGraphAssay`](CellGraphAssay-class.md) objects

- merge.data:

  Merge the data slots instead of just merging the counts (which
  requires renormalization); this is recommended if the same
  normalization approach was applied to all objects

- add.cell.ids:

  A character vector with sample names

- collapse:

  If TRUE, merge layers of the same name together

## Functions

- `RenameCells(CellGraphAssay)`: Rename cell IDs of a `CellGraphAssay`
  object

- `show(CellGraphAssay)`: Show method for `CellGraphAssay` objects

- `subset(CellGraphAssay)`: Subset a `CellGraphAssay` object

- `merge(CellGraphAssay)`: Merge two or more `CellGraphAssay` objects
  together
