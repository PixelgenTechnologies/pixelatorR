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
#' - 1x CD8 T cell (439026e80f97716f)
#' - 1x CD4 T cell (406b9e5d80941ca0)
#' - 1x NK cell (aae9e02d1dd14db5)
#' - 1x Monocyte (6add95ba33143eca)
#' - 1x B cell (3898b03349c6e28d)
#'
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
