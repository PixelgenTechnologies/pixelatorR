% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_data_mpx.R
\name{ReadMPX_counts}
\alias{ReadMPX_counts}
\title{Read a count matrix from a pxl file}
\usage{
ReadMPX_counts(filename, return_list = FALSE, verbose = TRUE)
}
\arguments{
\item{filename}{Path to a PXL file}

\item{return_list}{If TRUE, returns a list with the expression matrix and
the output from \code{h5read}}

\item{verbose}{Print messages}
}
\value{
A count matrix or a list if \code{return_list = TRUE}
}
\description{
Read a count matrix from a pxl file
}
\examples{
library(pixelatorR)

# Load example data
pxl_file <- minimal_mpx_pxl_file()
counts <- ReadMPX_counts(pxl_file)
counts[1:5, 1:5]

}
