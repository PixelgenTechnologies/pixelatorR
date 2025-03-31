#' Join Layers Together
#'
#' Join Layers Together
#'
#' @param object A \code{CellGraphAssay5} object.
#' @param layers A character vector of layer names to join.
#' @param new Name of new layers
#' @param ... Additional arguments passed to other methods
#'
#' @describeIn CellGraphAssay5-methods Join layers
#' @method JoinLayers CellGraphAssay5
#'
#' @return A \code{CellGraphAssay5} object with layers joined
#'
#' @examples
#' library(SeuratObject)
#' options(Seurat.object.assay.version = "v5")
#'
#' # Load example data as a Seurat object
#' pxl_file <- minimal_mpx_pxl_file()
#' seur_obj <- ReadMPX_Seurat(pxl_file)
#'
#' # Merge Seurat objects
#' seur_obj_merged <- merge(seur_obj, seur_obj)
#' cg_assay <- seur_obj_merged[["mpxCells"]]
#'
#' # The CellGraphAssay5 now has two count matrices
#' cg_assay
#'
#' # Join layers
#' cg_assay <- JoinLayers(cg_assay)
#'
#' # Now the CellGraphAssay5 has a single, merged count matrix
#' cg_assay
#'
#' # JoinLayers now also works on the Seurat object directly
#' seur_obj_merged <- JoinLayers(seur_obj_merged)
#' seur_obj_merged[["mpxCells"]]
#'
#' @export
#'
JoinLayers.CellGraphAssay5 <- function(
  object,
  layers = NULL,
  new = NULL,
  ...
) {
  assay5 <- as(object, "Assay5")
  assay5_joined <- JoinLayers(assay5, layers = layers, new = new, ... = ...)

  cg_assay5_joined <- as.CellGraphAssay5(
    x = assay5_joined,
    cellgraphs = CellGraphs(object),
    polarization = PolarizationScores(object),
    colocalization = ColocalizationScores(object),
    fs_map = FSMap(object)
  )

  return(cg_assay5_joined)
}


#' Join Layers Together
#'
#' Join Layers Together
#'
#' @param object A \code{PNAAssay5} object.
#' @param layers A character vector of layer names to join.
#' @param new Name of new layers
#' @param ... Additional arguments passed to other methods
#'
#' @describeIn PNAAssay5-methods Join layers
#' @method JoinLayers PNAAssay5
#'
#' @return A \code{PNAAssay5} object with layers joined
#'
#' @examples
#' library(SeuratObject)
#' options(Seurat.object.assay.version = "v5")
#' # Load example data as a Seurat object
#' pxl_file <- minimal_pna_pxl_file()
#' seur_obj <- ReadPNA_Seurat(pxl_file)
#'
#' # Merge Seurat objects
#' seur_obj_merged <- merge(seur_obj, seur_obj)
#' pna_assay <- seur_obj_merged[["PNA"]]
#'
#' # The PNAAssay5 now has two count matrices
#' pna_assay
#'
#' # Join layers
#' pna_assay <- JoinLayers(pna_assay)
#'
#' # Now the PNAAssay5 has a single, merged count matrix
#' pna_assay
#'
#' # JoinLayers now also works on the Seurat object directly
#' seur_obj_merged <- JoinLayers(seur_obj_merged)
#' seur_obj_merged[["PNA"]]
#'
#' @export
#'
JoinLayers.PNAAssay5 <- function(
  object,
  layers = NULL,
  new = NULL,
  ...
) {
  assay5 <- as(object, "Assay5")
  assay5_joined <- JoinLayers(assay5, layers = layers, new = new, ... = ...)

  cg_assay5_joined <- as.PNAAssay5(
    x = assay5_joined,
    cellgraphs = CellGraphs(object),
    proximity = ProximityScores(object),
    fs_map = FSMap(object)
  )

  return(cg_assay5_joined)
}
