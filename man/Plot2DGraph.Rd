% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph_layout_visualization.R
\name{Plot2DGraph}
\alias{Plot2DGraph}
\title{Plot 2D graph layouts}
\usage{
Plot2DGraph(
  object,
  cells = NULL,
  marker = NULL,
  assay = NULL,
  layout_method = c("pmds", "wpmds", "fr", "kk", "drl"),
  colors = c("lightgrey", "mistyrose", "red", "darkred"),
  map_nodes = TRUE,
  map_edges = FALSE,
  node_size = 0.5,
  edge_width = 0.3,
  ...
)
}
\arguments{
\item{object}{A \code{Seurat} object}

\item{cells}{A character vector with cell IDs}

\item{marker}{Name of a marker to colors nodes/edges by}

\item{assay}{Name of assay to pull data from}

\item{layout_method}{Select appropriate layout previously computed with
\code{\link{ComputeLayout}}}

\item{colors}{A character vector of colors to color marker counts by}

\item{map_nodes, map_edges}{Should nodes and/or edges be mapped?}

\item{node_size}{Size of nodes}

\item{edge_width}{Not yet implemented TODO}

\item{...}{Not yet implemented}
}
\value{
An object of class \code{patchwork}
}
\description{
Plot 2D component graph layouts computed with \code{\link{ComputeLayout}} and
optionally color nodes by certain attributes. Edges can also be visualized
by setting \code{map_edges}; however, since component graphs tend to be very
large, this can take a long time to draw.
}
\examples{
library(pixelatorR)
pxl_file <- system.file("extdata/PBMC_10_cells",
                        "Sample01_test.pxl",
                        package = "pixelatorR")

seur <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)
seur <- LoadCellGraphs(seur, load_as = "Anode", cells = colnames(seur)[1:10])
seur[["mpxCells"]] <- KeepLargestComponent(seur[["mpxCells"]])
seur <- ComputeLayout(seur, layout_method = "pmds", dim = 2)

Plot2DGraph(seur, cells = colnames(seur)[1], marker = "HLA-ABC")

}