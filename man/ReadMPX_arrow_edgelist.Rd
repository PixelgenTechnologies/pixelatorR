% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_edgelist.R
\name{ReadMPX_arrow_edgelist}
\alias{ReadMPX_arrow_edgelist}
\title{Read edgelists from a PXL file containing MPX data}
\usage{
ReadMPX_arrow_edgelist(pxl_file, edge_list_file = NULL, verbose = TRUE, ...)
}
\arguments{
\item{pxl_file}{Path to a PXL file}

\item{edge_list_file}{Path to the output edgelist.parquet file}

\item{verbose}{Print messages}

\item{...}{Parameters passed to other methods}
}
\value{
Nothing. The edgelist is saved to a parquet file
set with \code{edge_list_file}
}
\description{
This function uses arrow to read edgelists from one or several
PXL files. The edgelists are stored in parquet files in
the \code{outdir} directory which can be modified on disk.
}
\examples{
library(pixelatorR)

# Load example data
pxl_file <- minimal_mpx_pxl_file()
edgelist_arrow <- ReadMPX_arrow_edgelist(pxl_file)
edgelist_arrow

}
