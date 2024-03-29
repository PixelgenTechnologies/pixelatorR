% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_and_write_parquet.R
\name{ReadMPX_arrow_edgelist}
\alias{ReadMPX_arrow_edgelist}
\title{Read edgelists from .pxl files}
\usage{
ReadMPX_arrow_edgelist(
  path,
  outdir = NULL,
  return_list = FALSE,
  overwrite = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{path}{Path to a .pxl file or a directory containing .parquet files}

\item{outdir}{A path to an existing directory}

\item{return_list}{If TRUE, returns a list with an \code{ArrowObject},
a vector with sample IDs and a the path to the parquet output directory}

\item{overwrite}{If the output directory already exists, set this parameter
to TRUE to overwrite the directory}

\item{verbose}{Print messages}
}
\value{
An \code{ArrowObject}
}
\description{
This function uses arrow to read edgelists from one or several
.pxl files. The edgelists are stored in parquet files in
the \code{outdir} directory which can be modified on disk.
}
\examples{

library(pixelatorR)
# Set arrow data output directory to temp for tests
options(pixelatorR.arrow_outdir = tempdir())

# Load example data
pxl_file <- system.file("extdata/five_cells",
                        "five_cells.pxl",
                        package = "pixelatorR")
edgelist_arrow <- ReadMPX_arrow_edgelist(pxl_file, overwrite = TRUE)
edgelist_arrow

}
