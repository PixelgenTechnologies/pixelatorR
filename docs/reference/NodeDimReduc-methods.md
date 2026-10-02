# NodeDimReduc Methods

Methods for [`NodeDimReduc`](NodeDimReduc-class.md) objects

## Usage

``` r
# S4 method for class 'NodeDimReduc'
show(object)

# S3 method for class 'NodeDimReduc'
Embeddings(object, ...)

# S3 method for class 'NodeDimReduc'
Loadings(object, projected = FALSE, ...)

# S3 method for class 'NodeDimReduc'
Stdev(object, ...)

# S3 method for class 'NodeDimReduc'
Key(object, ...)

# S3 method for class 'NodeDimReduc'
Cells(x, ...)
```

## Arguments

- object:

  A [`NodeDimReduc`](NodeDimReduc-class.md) object

- ...:

  Currently not used

- projected:

  Ignored; included for compatibility with
  [`Loadings`](https://satijalab.github.io/seurat-object/reference/Loadings.html)

- x:

  A [`NodeDimReduc`](NodeDimReduc-class.md) object

## Functions

- `show(NodeDimReduc)`: Show a `NodeDimReduc` object
