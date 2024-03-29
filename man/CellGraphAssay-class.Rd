% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CellGraphAssay.R
\docType{class}
\name{CellGraphAssay-class}
\alias{CellGraphAssay-class}
\alias{CellGraphAssay}
\title{The CellGraphAssay class}
\description{
The CellGraphAssay object is an extended \code{\link[SeuratObject]{Assay}}
for the storage and analysis of mpx single-cell data.
}
\section{Slots}{

\describe{
\item{\code{cellgraphs}}{A named list of \code{\link{CellGraph}} objects}

\item{\code{polarization}}{A \code{tbl_df} with polarization scores}

\item{\code{colocalization}}{A \code{tbl_df} with colocalization scores}

\item{\code{arrow_dir}}{A character giving the name of the directory where the edgelist
parquet file(s) are stored}

\item{\code{arrow_data}}{An \code{R6} class object with an arrow Dataset}
}}

\concept{assay}
