% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/graph_layout.R
\name{ComputeLayout}
\alias{ComputeLayout}
\alias{ComputeLayout.tbl_graph}
\alias{ComputeLayout.CellGraph}
\alias{ComputeLayout.CellGraphAssay}
\alias{ComputeLayout.Seurat}
\title{Compute a graph layout}
\usage{
ComputeLayout(object, ...)

\method{ComputeLayout}{tbl_graph}(
  object,
  layout_method = c("pmds", "fr", "kk", "drl"),
  dim = 2,
  normalize_layout = FALSE,
  project_on_unit_sphere = FALSE,
  k = 0,
  pivots = 100,
  seed = 123,
  custom_layout_function = NULL,
  custom_layout_function_args = NULL,
  ...
)

\method{ComputeLayout}{CellGraph}(
  object,
  layout_method = c("pmds", "fr", "kk", "drl"),
  dim = 2,
  normalize_layout = FALSE,
  project_on_unit_sphere = FALSE,
  k = 0,
  pivots = 100,
  seed = 123,
  custom_layout_function = NULL,
  custom_layout_function_args = NULL,
  custom_layout_name = "custom",
  ...
)

\method{ComputeLayout}{CellGraphAssay}(
  object,
  layout_method = c("pmds", "fr", "kk", "drl"),
  dim = 2,
  normalize_layout = FALSE,
  project_on_unit_sphere = FALSE,
  k = 0,
  pivots = 100,
  seed = 123,
  verbose = TRUE,
  custom_layout_function = NULL,
  custom_layout_function_args = NULL,
  custom_layout_name = "custom",
  ...
)

\method{ComputeLayout}{Seurat}(
  object,
  assay = NULL,
  layout_method = c("pmds", "fr", "kk", "drl"),
  dim = 2,
  normalize_layout = FALSE,
  project_on_unit_sphere = FALSE,
  k = 0,
  pivots = 100,
  seed = 123,
  verbose = TRUE,
  custom_layout_function = NULL,
  custom_layout_function_args = NULL,
  custom_layout_name = "custom",
  ...
)
}
\arguments{
\item{object}{An object}

\item{...}{Additional parameters passed to other methods}

\item{layout_method}{The method for calculating the graph layout; PMDS,
Fruchterman-Reingold (fr), Kamada-Kawai (kk), drl.}

\item{dim}{An integer specifying the dimensions of the layout (2 or 3)}

\item{normalize_layout}{Logical specifying whether the coordinate system
should be centered at origo and the coordinates scaled such that their median
length (euclidean norm) is 1.}

\item{project_on_unit_sphere}{Should the resulting layout be projected onto
a unit sphere?}

\item{k}{The size of the neighborhood from which to pool counts from in
the UPIA antibody count table. 0 is recommended.}

\item{pivots}{Only used for "pmds" graph layout. See \code{?layout_with_pmds}
for details}

\item{seed}{Set seed for reproducibility}

\item{custom_layout_function}{A custom function for layout computation. The
function should take a \code{tbl_graph} object and a \code{dim} value as
input and return a matrix of dimensions NxD, where N is the number of nodes
and D is equal to \code{dim}. Note that this will override the \code{layout_method}.}

\item{custom_layout_function_args}{A list of arguments passed to \code{custom_layout_function}.
The \code{dim} is automatically passed to \code{custom_layout_function} and should not be
included in \code{custom_layout_function_args}.}

\item{custom_layout_name}{A name for the layout computed with
\code{custom_layout_function}. Should not be one of "pmds", "fr",
"kk" or "drl".}

\item{verbose}{Print messages}

\item{assay}{Name of assay to compute layouts for}
}
\value{
An object containing a graph layout
}
\description{
Calculates a graph layout for a component's edgelist, and outputs a list
with the bipartite graph, layout, and antibody counts per A pixel.
}
\examples{
library(pixelatorR)
library(dplyr)
# Set arrow data output directory to temp for tests
options(pixelatorR.arrow_outdir = tempdir())

pxl_file <- system.file("extdata/five_cells",
                        "five_cells.pxl",
                        package = "pixelatorR")

# Load example data
seur <- ReadMPX_Seurat(pxl_file, return_cellgraphassay = TRUE, overwrite = TRUE)

# Load 1 cellgraph
seur <- LoadCellGraphs(seur, cells = colnames(seur)[1],
                       load_as = "Anode", force = TRUE)

# Get CellGraph
cg <- CellGraphs(seur)[[colnames(seur)[1]]]

# Get tbl_graph object
tbl_graph <- slot(cg, name = "cellgraph")

# Compute layout for a tbl_graph
layout <- ComputeLayout(tbl_graph, layout_method = "fr")
layout \%>\% head()


# Compute layout for a CellGraph
cg <- ComputeLayout(cg, layout_method = "fr")


# Compute layout for a CellGraphAssay
cg_assay <- ComputeLayout(seur[["mpxCells"]], layout_method = "fr")


# Compute layout for a Seurat object
seur <- ComputeLayout(seur, layout_method = "fr")

}
\seealso{
\code{\link[=center_layout_coordinates]{center_layout_coordinates()}} for centering layout coordinates,
\code{\link[=normalize_layout_coordinates]{normalize_layout_coordinates()}} for normalizing layout coordinates and
\code{\link[=project_layout_coordinates_on_unit_sphere]{project_layout_coordinates_on_unit_sphere()}} for projecting layout coordinates onto a unit sphere.
}
