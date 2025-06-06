% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph_layout_visualization.R
\name{Plot3DGraph}
\alias{Plot3DGraph}
\title{Plot 3D graph layouts}
\usage{
Plot3DGraph(
  object,
  cell_id,
  marker = NULL,
  assay = NULL,
  layout_method = c("wpmds_3d", "pmds_3d"),
  project = FALSE,
  aspectmode = c("data", "cube"),
  colors = c("lightgrey", "mistyrose", "red", "darkred"),
  showgrid = TRUE,
  log_scale = TRUE,
  node_size = 2,
  show_Bnodes = FALSE,
  ...
)
}
\arguments{
\item{object}{A \code{Seurat} object}

\item{cell_id}{ID of component to visualize}

\item{marker}{Name of marker to color the nodes by}

\item{assay}{Name of assay to pull data from}

\item{layout_method}{Select appropriate layout previously computed with
\code{\link{ComputeLayout}}}

\item{project}{Project the nodes onto a sphere. Default FALSE}

\item{aspectmode}{Set aspect ratio to one of "data" or "cube".
If "cube", this scene's axes are drawn as a cube, regardless of the axes' ranges.
If "data", this scene's axes are drawn in proportion with the axes' ranges.

Default "data"}

\item{colors}{Color the nodes expressing a marker. Must be a character vector
with at least two color names.}

\item{showgrid}{Show the grid lines. Default TRUE}

\item{log_scale}{Convert node counts to log-scale with \code{logp}}

\item{node_size}{Size of nodes}

\item{show_Bnodes}{Should B nodes be included in the visualization?
This option is only applicable to bipartite graphs.}

\item{...}{Additional parameters passed to \code{plot_ly}}
}
\value{
A interactive 3D plot of a component graph layout as a \code{plotly} object
}
\description{
Plot a 3D component graph layout computed with \code{\link{ComputeLayout}} and
color nodes by a marker.
}
\examples{
library(pixelatorR)

# MPX
pxl_file <- minimal_mpx_pxl_file()
seur <- ReadMPX_Seurat(pxl_file)
seur <- LoadCellGraphs(seur, cells = colnames(seur)[5])
seur <- ComputeLayout(seur, layout_method = "wpmds", dim = 3, pivots = 50)
Plot3DGraph(seur, cell_id = colnames(seur)[5], marker = "CD50", layout_method = "wpmds_3d")

# PNA
pxl_file <- minimal_pna_pxl_file()
seur <- ReadPNA_Seurat(pxl_file)
seur <- LoadCellGraphs(seur, cells = colnames(seur)[1], add_layouts = TRUE)
Plot3DGraph(seur, cell_id = colnames(seur)[1], marker = "CD16", layout_method = "wpmds_3d")

}
