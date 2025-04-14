#' Five Cells MPX Test data
#'
#' MPX data set of 5 immune cells (link to [data resource](https://software.pixelgen.com/datasets))
#' \itemize{
#'    \item{1x CD3 capping stimulated CD4 T cell (RCVCMP0000217) \cr\cr
#'    \strong{Data set:} CD3 Stimulated Human PBMCs SCSP v1.0 Kit, Immunology Panel I}
#'    \item{1x B cell (RCVCMP0000118)\cr
#'    1x CD4 T cell (RCVCMP0000487)\cr
#'    1x CD8 T cell (RCVCMP0000655) \cr\cr
#'    \strong{Data set:} 1k Human PBMCs, SCSP v1.0 Kit, Immunology Panel I}
#'    \item{1x Migratory CD8 T cell with uropod (RCVCMP0000263) \cr\cr
#'    \strong{Data set:} Uropod formation on T cells, SCSP v1.0 Kit, Immunology Panel I}
#' }
#'
#' @name mpx_dataset
#' @family datasets
NULL

#' Five Cells PNA Test data
#'
#' PNA data set of 5 PBMC immune cells
#'
#' - 1x CD16+ Mono cell (0a45497c6bfbfb22)
#' - 1x pDC cell (2708240b908e2eba)
#' - 3x CD4 T cell (c3c393e9a17c1981, d4074c845bb62800, efe0ed189cb499fc)
#'
#' @name pna_dataset
#' @family datasets
NULL

#' Get path to test MPX PXL file
#'
#' @rdname mpx_dataset
#'
#' @return Path to the PXL file containing a minimal MPX data set used in tests
#'
#' @export
#'
minimal_mpx_pxl_file <- function() {
  system.file("extdata/five_cells",
    "five_cells.pxl",
    package = "pixelatorR"
  )
}

#' Get path to test PNA PXL file
#'
#' @rdname pna_dataset
#'
#' @return Path to the PXL file containing a minimal PNA data set used in tests
#'
#' @export
#'
minimal_pna_pxl_file <- function() {
  system.file("extdata/five_cells",
    "minimal_PNA_PBMC.pxl",
    package = "pixelatorR"
  )
}
