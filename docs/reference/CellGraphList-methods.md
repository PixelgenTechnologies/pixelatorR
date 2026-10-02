# CellGraphList Methods

Methods for [`CellGraphList`](CellGraphList.md) objects. Subsetting,
concatenation, and replacement type-check elements against
[`CellGraph`](CellGraph-class.md). Unloaded graphs may be stored as
`NULL`; `x[[i]] <- NULL` keeps the name and stores `NULL` rather than
dropping the element.

## Usage

``` r
# S3 method for class 'CellGraphList'
print(x, ...)

# S3 method for class 'CellGraphList'
x[i, ...]

# S3 method for class 'CellGraphList'
x[i] <- value

# S3 method for class 'CellGraphList'
x[[i]] <- value

# S3 method for class 'CellGraphList'
names(x) <- value

# S3 method for class 'CellGraphList'
c(...)

# S3 method for class 'CellGraphList'
as.list(x, ...)
```

## Arguments

- x:

  A [`CellGraphList`](CellGraphList.md) object

- ...:

  Currently not used

- i:

  Index to extract or replace

- value:

  A [`CellGraph`](CellGraph-class.md), `NULL`, or a list of those

## Value

`[`, `[<-`, `[[<-`, `names<-`, and `c`: a `CellGraphList`. `as.list`: a
named list. `print`: `x`, invisibly.

## Functions

- `print(CellGraphList)`: Print a `CellGraphList`

- `[`: Subset a `CellGraphList`. Unknown character names raise an error.

- `` `[`(CellGraphList) <- value ``: Replace a subset of graphs. `NULL`
  unloads the selected cells without dropping their names.

- `` `[[`(CellGraphList) <- value ``: Replace a single graph. `NULL`
  unloads that cell without dropping its name.

- `names(CellGraphList) <- value`: Set names. Names must be unique and
  non-missing.

- `c(CellGraphList)`: Concatenate `CellGraphList` objects with lists of
  `CellGraph` or `NULL`. A bare `CellGraph` uses the argument name, or
  `CellGraph1`, `CellGraph2`, ... when unnamed.

- `as.list(CellGraphList)`: Convert to a named list

## See also

[`CellGraphList`](CellGraphList.md)

## Examples

``` r
library(pixelatorR)
library(tidygraph)
library(dplyr)

edges <- tibble(from = c("a", "b"), to = c("b", "c"))
g <- as_tbl_graph(edges, directed = FALSE) %N>%
  mutate(node_type = c("umi1", "umi2", "umi1"))
attr(g, "type") <- "bipartite"
cg <- CreateCellGraphObject(cellgraph = g)
cgl <- CreateCellGraphList(list(cell_1 = cg, cell_2 = cg))

# Print and subset
print(cgl)
#> A CellGraphList with 2 loaded CellGraph object(s) out of 2 
#> Names: cell_1, cell_2  
cgl[1]
#> A CellGraphList with 1 loaded CellGraph object(s) out of 1 
#> Names: cell_1  

# Unload a graph without dropping its name
cgl[[1]] <- NULL

# Concatenate while keeping the class
c(cgl[1], cgl[2])
#> A CellGraphList with 1 loaded CellGraph object(s) out of 2 
#> Names: cell_1, cell_2  
```
