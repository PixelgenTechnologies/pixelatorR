# Summarise a cell plot recipe

Creates a compact summary of a recipe without building or drawing it.

## Usage

``` r
# S3 method for class 'cell_plot'
summary(object, ...)
```

## Arguments

- object:

  A `cell_plot` recipe.

- ...:

  Additional arguments. Currently not used.

## Value

A `summary.cell_plot` object. Printing writes one line per field. The
underlying value is a character vector describing the recipe.
