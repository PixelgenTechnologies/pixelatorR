% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PNAAssay.R
\docType{class}
\name{PNAAssay5-class}
\alias{PNAAssay5-class}
\alias{PNAAssay5}
\title{The PNAAssay5 class}
\description{
The PNAAssay5 object is an extended \code{\link[SeuratObject]{Assay5}}
for storing PNA single-cell data.
}
\details{
Compared to the \code{\link[SeuratObject]{Assay5}} class, the PNAAssay5
class has three additional slots:
\itemize{
\item \code{cellgraphs}: A named list of \code{\link{CellGraph}} objects
\item \code{proximity}: A \code{tbl_df} with proximity scores
\item \code{fs_map}: A \code{tbl_df} with information on source PXL file
}
}
\section{Slots}{

\describe{
\item{\code{cellgraphs}}{A named list of \code{\link{CellGraph}} objects}

\item{\code{proximity}}{A \code{tbl_df} with proximity scores}

\item{\code{fs_map}}{A \code{tbl_df} with information on source PXL file
paths, sample IDs, and component IDs}
}}

\concept{assay}
