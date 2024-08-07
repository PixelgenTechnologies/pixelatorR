% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/CellGraphAssay.R, R/objects.R
\name{CellGraphs}
\alias{CellGraphs}
\alias{CellGraphs<-}
\alias{CellGraphs.MPXAssay}
\alias{CellGraphs<-.MPXAssay}
\alias{CellGraphs.Seurat}
\alias{CellGraphs<-.Seurat}
\title{CellGraphs}
\usage{
CellGraphs(object, ...)

CellGraphs(object, ...) <- value

\method{CellGraphs}{MPXAssay}(object, ...)

\method{CellGraphs}{MPXAssay}(object, ...) <- value

\method{CellGraphs}{Seurat}(object, ...)

\method{CellGraphs}{Seurat}(object, ...) <- value
}
\arguments{
\item{object}{A \code{Seurat}, \code{CellGraphAssay} or \code{CellGraphAssay} object}

\item{...}{Arguments passed to other methods}

\item{value}{A list of \code{\link{CellGraph}} objects}
}
\value{
Returns a list of \code{\link{CellGraph}} objects. If there are
no \code{\link{CellGraph}} objects present, returns an empty named list.
}
\description{
Get and set \code{\link{CellGraph}} lists for different objects.
}
\examples{

library(pixelatorR)
library(dplyr)
library(tidygraph)

pxl_file <- system.file("extdata/five_cells",
  "five_cells.pxl",
  package = "pixelatorR"
)
counts <- ReadMPX_counts(pxl_file)
edgelist <- ReadMPX_item(pxl_file, items = "edgelist")
components <- colnames(counts)
edgelist_split <-
  edgelist \%>\%
  select(upia, upib, component) \%>\%
  distinct() \%>\%
  group_by(component) \%>\%
  group_split() \%>\%
  setNames(nm = components)

# Convert data into a list of CellGraph objects
bipartite_graphs <- lapply(edgelist_split, function(x) {
  x <- x \%>\% as_tbl_graph(directed = FALSE)
  x <- x \%>\% mutate(node_type = case_when(name \%in\% edgelist$upia ~ "A", TRUE ~ "B"))
  attr(x, "type") <- "bipartite"
  CreateCellGraphObject(cellgraph = x)
})

# CellGraphs getter CellGraphAssay
# ---------------------------------

# Create CellGraphAssay
cg_assay <- CreateCellGraphAssay(counts = counts, cellgraphs = bipartite_graphs)
cg_assay

# Get cellgraphs from a CellGraphAssay object
CellGraphs(cg_assay)


# CellGraphs setter CellGraphAssay
# ---------------------------------

# Set cellgraphs in a CellGraphAssay object
CellGraphs(cg_assay) <- cg_assay@cellgraphs


# CellGraphs getter Seurat
# ---------------------------------
se <- ReadMPX_Seurat(pxl_file)

# Get cellgraphs from a Seurat object
CellGraphs(se)

# CellGraphs setter Seurat
# ---------------------------------

# Set cellgraphs in a Seurat object
CellGraphs(se) <- cg_assay@cellgraphs

}
\seealso{
\code{\link[=PolarizationScores]{PolarizationScores()}} and \code{\link[=ColocalizationScores]{ColocalizationScores()}} for getting/setting spatial metrics
}
\concept{assay}
\concept{cellgraphs}
