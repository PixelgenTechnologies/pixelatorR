% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CellGraph.R
\name{CreateCellGraphObject}
\alias{CreateCellGraphObject}
\title{Create a CellGraph object}
\usage{
CreateCellGraphObject(cellgraph, counts = NULL, layout = NULL, verbose = FALSE)
}
\arguments{
\item{cellgraph}{A \code{tbl_graph} object representing an mpx single-cell graph}

\item{counts}{A \code{dgCMatrix} with marker counts}

\item{layout}{A \code{tbl_df} object with cell layout(s)}

\item{verbose}{Print messages}
}
\value{
A \code{CellGraph} object
}
\description{
Create a CellGraph object
}
\examples{

library(pixelatorR)
library(dplyr)
library(tidygraph)

edge_list <-
  ReadMPX_item(
    minimal_mpx_pxl_file(),
    items = "edgelist"
  )
bipart_graph <-
  edge_list \%>\%
  select(upia, upib, marker) \%>\%
  distinct() \%>\%
  as_tbl_graph(directed = FALSE) \%>\%
  mutate(node_type = case_when(name \%in\% edge_list$upia ~ "A", TRUE ~ "B"))
attr(bipart_graph, "type") <- "bipartite"

cg <- CreateCellGraphObject(cellgraph = bipart_graph)
cg

}
\concept{cellgraph}
