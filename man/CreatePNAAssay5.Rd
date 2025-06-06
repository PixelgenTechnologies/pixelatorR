% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PNAAssay.R
\name{CreatePNAAssay5}
\alias{CreatePNAAssay5}
\title{Create a PNAAssay5 object}
\usage{
CreatePNAAssay5(
  counts,
  cellgraphs,
  proximity = NULL,
  fs_map = NULL,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{counts}{Unnormalized data (raw counts)}

\item{cellgraphs}{A named list with \code{\link{CellGraph}} objects}

\item{proximity}{A \code{tbl_df} with proximity scores}

\item{fs_map}{A \code{tbl_df} with information on source PXL file
paths, sample IDs, and component IDs}

\item{verbose}{Print messages}

\item{...}{Additional arguments passed to \code{\link[SeuratObject]{CreateAssay5Object}}}
}
\value{
A \code{PNAAssay5} object
}
\description{
Create a \code{\link{PNAAssay5}} object from a count matrix and a list
of \code{\link{CellGraph}} objects. The expected format of the input
matrix is features x cells. Optionally, a \code{tbl_df} with proximity
scores and a \code{tbl_df} with information on source PXL file paths
can be provided.
}
\examples{
library(pixelatorR)
library(dplyr)

pxl_file <- minimal_pna_pxl_file()
counts <- ReadPNA_counts(pxl_file)
pna_assay5 <- CreatePNAAssay5(
  counts = counts,
  cellgraphs = rep(list(NULL), ncol(counts)) \%>\%
    setNames(colnames(counts))
)
pna_assay5

}
\concept{assay}
