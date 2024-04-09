#' @include generics.R
NULL

#' @rdname PolarizationScoresToAssay
#' @method PolarizationScoresToAssay data.frame
#'
#' @examples
#' library(pixelatorR)
#' library(SeuratObject)
#' # Set arrow data output directory to temp for tests
#' options(pixelatorR.arrow_outdir = tempdir())
#'
#' # Load example data as a Seurat object
#' pxl_file <- system.file("extdata/five_cells",
#'   "five_cells.pxl",
#'   package = "pixelatorR"
#' )
#' pol_scores <- ReadMPX_polarization(pxl_file)
#'
#' # PolarizationScoresToAssay returns a matrix for a tbl_df
#' pol_scores_mat <- PolarizationScoresToAssay(pol_scores)
#' pol_scores_mat[1:4, 1:5]
#'
#' @export
#'
PolarizationScoresToAssay.data.frame <- function (
  object,
  values_from = c("morans_z", "morans_i"),
  ...
) {

  # Validate input
  values_from <- match.arg(values_from, choices = c("morans_z", "morans_i"))
  stopifnot(
    "'component' must be present in input table" = "component" %in% colnames(object),
    "'marker' must be present in input table" = "marker" %in% colnames(object),
    "'values_from' must be present in input table" = values_from %in% colnames(object)
  )

  # Cast data.frame to wide format
  pol_scores_wide_format <- object %>%
    pivot_wider(id_cols = "marker",
                names_from = "component",
                values_from = all_of(values_from),
                values_fill = 0) %>%
    data.frame(row.names = 1) %>%
    as.matrix()

  # Replace missing values with 0
  pol_scores_wide_format[is.na(pol_scores_wide_format)] <- 0

  return(pol_scores_wide_format)
}


#' @rdname PolarizationScoresToAssay
#' @method PolarizationScoresToAssay CellGraphAssay
#'
#' @examples
#' # Create a Seurat object
#' seur <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)
#'
#' # Fetch CellGraphAssay and create new polarization
#' # scores Assay
#' cg_assay <- seur[["mpxCells"]]
#' class(cg_assay)
#' pol_assay <- PolarizationScoresToAssay(cg_assay)
#' class(pol_assay)
#'
#' @export
#'
PolarizationScoresToAssay.CellGraphAssay <- function (
  object,
  values_from = c("morans_z", "morans_i"),
  ...
) {

  pol_matrix <- .create_spatial_metric_matrix(object,
                                              values_from = values_from,
                                              metric = "polarization")

  # Create Assay from filled matrix
  assay <- CreateAssayObject(data = pol_matrix)

  return(assay)
}

#' @rdname PolarizationScoresToAssay
#' @method PolarizationScoresToAssay CellGraphAssay5
#'
#' @export
#'
PolarizationScoresToAssay.CellGraphAssay5 <- function (
  object,
  values_from = c("morans_z", "morans_i"),
  ...
) {

  pol_matrix <- .create_spatial_metric_matrix(object,
                                              values_from = values_from,
                                              metric = "polarization")

  # Create Assay from filled matrix
  assay <- CreateAssay5Object(data = pol_matrix)

  return(assay)
}


#' @param assay Name of the \code{\link{CellGraphAssay}} to pull polarization scores from
#' @param new_assay Name of the \code{Assay} to store the polarization scores in
#'
#' @rdname PolarizationScoresToAssay
#' @method PolarizationScoresToAssay Seurat
#'
#' @examples
#' # Convert polarization scores within a Seurat object
#' seur <- PolarizationScoresToAssay(seur)
#'
#' # After conversion, we now have a new assay called "polarization"
#' seur[["polarization"]]
#'
#' # Switch default assay to polarization
#' DefaultAssay(seur) <- "polarization"
#'
#' # Visualize polarization scores with Seurat
#' # VlnPlot(seur, features = "CD19") +
#' #   labs(y = "Polarization score")
#'
#' @export
#'
PolarizationScoresToAssay.Seurat <- function (
  object,
  assay = NULL,
  new_assay = NULL,
  values_from = c("morans_z", "morans_i"),
  ...
) {

  # Use default assay if assay = NULL
  if (!is.null(assay)) {
    stopifnot(
      "'assay' must be a character of length 1" =
        is.character(assay) &&
        (length(assay) == 1)
    )
  } else {
    # Use default assay if assay = NULL
    assay <- DefaultAssay(object)
  }

  # Use "polarization" if new_assay = NULL
  if (!is.null(new_assay)) {
    stopifnot(
      "'new_assay' must be a character of length 1" =
        is.character(new_assay) &&
        (length(new_assay) == 1)
    )
  } else {
    # Use default assay if assay = NULL
    new_assay <- "polarization"
  }

  # fetch CellGraphAssay
  cg_assay <- object[[assay]]

  # Validate input
  values_from <- match.arg(values_from, choices = c("morans_z", "morans_i"))

  # Create new assay
  pol_assay <- PolarizationScoresToAssay(cg_assay, values_from, ...)

  # Place new assay in Seurat object
  object[[new_assay]] <- pol_assay

  return(object)
}


#' @rdname ColocalizationScoresToAssay
#' @method ColocalizationScoresToAssay data.frame
#'
#' @examples
#' library(pixelatorR)
#' library(SeuratObject)
#' # Set arrow data output directory to temp for tests
#' options(pixelatorR.arrow_outdir = tempdir())
#'
#' # Load example data as a Seurat object
#' pxl_file <- system.file("extdata/five_cells",
#'   "five_cells.pxl",
#'   package = "pixelatorR"
#' )
#' col_scores <- ReadMPX_colocalization(pxl_file)
#'
#' # ColocalizationScoresToAssay returns a matrix for a tbl_df
#' col_scores_mat <- ColocalizationScoresToAssay(col_scores)
#' col_scores_mat[1:4, 1:5]
#'
#' @export
#'
ColocalizationScoresToAssay.data.frame <- function (
  object,
  values_from = c("pearson_z", "pearson"),
  ...
) {

  # Validate input
  values_from <- match.arg(values_from,
                           choices = c("pearson_z", "pearson"))
  stopifnot(
    "'component' must be present in input table" =
      "component" %in% colnames(object),
    "'marker_1' must be present in input table" =
      "marker_1" %in% colnames(object),
    "'marker_2' must be present in input table" =
      "marker_2" %in% colnames(object),
    "'values_from' must be present in input table" =
      values_from %in% colnames(object)
  )

  # Cast data.frame to wide format
  col_scores_wide_format <- object %>%
    pivot_wider(id_cols = c("marker_1", "marker_2"),
                names_from = "component",
                values_from = all_of(values_from),
                values_fill = 0) %>%
    dplyr::filter(marker_1 != marker_2) %>%
    unite(marker_1, marker_2, col = "pair", sep = "-") %>%
    data.frame(row.names = 1) %>%
    as.matrix()

  # Replace missing values with 0
  col_scores_wide_format[is.na(col_scores_wide_format)] <- 0

  return(col_scores_wide_format)
}


#' @rdname ColocalizationScoresToAssay
#' @method ColocalizationScoresToAssay CellGraphAssay
#'
#' @examples
#' # Create a Seurat object
#' seur <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)
#'
#' # Fetch CellGraphAssay and create new polarization
#' # scores Assay
#' cg_assay <- seur[["mpxCells"]]
#' class(cg_assay)
#' col_assay <- ColocalizationScoresToAssay(cg_assay)
#' class(col_assay)
#'
#' @export
#'
ColocalizationScoresToAssay.CellGraphAssay <- function (
  object,
  values_from = c("pearson_z", "pearson"),
  ...
) {

  coloc_matrix <- .create_spatial_metric_matrix(object,
                                                values_from = values_from,
                                                metric = "colocalization")

  # Create Assay from filled matrix
  assay <- CreateAssayObject(data = coloc_matrix)

  return(assay)
}


#' @rdname ColocalizationScoresToAssay
#' @method ColocalizationScoresToAssay CellGraphAssay5
#'
#' @export
#'
ColocalizationScoresToAssay.CellGraphAssay5 <- function (
  object,
  values_from = c("pearson_z", "pearson"),
  ...
) {

  coloc_matrix <- .create_spatial_metric_matrix(object,
                                                values_from = values_from,
                                                metric = "colocalization")

  # Create Assay from filled matrix
  assay <- CreateAssay5Object(data = coloc_matrix)

  return(assay)
}


#' @param assay Name of the \code{\link{CellGraphAssay}} to pull polarization scores from
#' @param new_assay Name of the \code{Assay} to store the polarization scores in
#'
#' @rdname ColocalizationScoresToAssay
#' @method ColocalizationScoresToAssay Seurat
#'
#' @examples
#' # Convert colocalization scores within a Seurat object
#' seur <- ColocalizationScoresToAssay(seur)
#'
#' # After conversion, we now have a new assay called "colocalization"
#' seur[["colocalization"]]
#'
#' # Switch default assay to polarization
#' DefaultAssay(seur) <- "colocalization"
#'
#' # Visualize colocalization scores with Seurat
#' # VlnPlot(seur, features = "CD19") +
#' #   ggplot2::labs(y = "Colocalization score")
#'
#' @export
#'
ColocalizationScoresToAssay.Seurat <- function (
  object,
  assay = NULL,
  new_assay = NULL,
  values_from = c("pearson_z", "pearson"),
  ...
) {

  # Use default assay if assay = NULL
  if (!is.null(assay)) {
    stopifnot(
      "'assay' must be a character of length 1" =
        is.character(assay) &&
        (length(assay) == 1)
    )
  } else {
    # Use default assay if assay = NULL
    assay <- DefaultAssay(object)
  }

  # Use "polarization" if new_assay = NULL
  if (!is.null(new_assay)) {
    stopifnot(
      "'new_assay' must be a character of length 1" =
        is.character(new_assay) &&
        (length(new_assay) == 1)
    )
  } else {
    # Use default assay if assay = NULL
    new_assay <- "colocalization"
  }

  # fetch CellGraphAssay
  cg_assay <- object[[assay]]

  # Validate input
  values_from <- match.arg(values_from,
                           choices = c("pearson_z", "pearson"))

  # Create new assay
  col_assay <- ColocalizationScoresToAssay(cg_assay, values_from, ...)

  # Place new assay in Seurat object
  object[[new_assay]] <- col_assay

  return(object)
}

#' Utility function to fetch spatial metrics from
#' a CellGraphAssay(5) object and convert them to a matrix
#'
#' @noRd
#'
.create_spatial_metric_matrix <- function (
  object,
  values_from,
  metric,
  ...
) {
  # fetch polarization scores
  spatial_metric_table <- slot(object, name = metric)
  if (is.null(spatial_metric_table)) {
    abort(glue("'{metric}' scores are missing from '{class(object)}'"))
  }

  # Validate input
  values_from <- match.arg(values_from, choices = switch(
    metric,
    "polarization" = c("morans_z", "morans_i"),
    "colocalization" = c("pearson_z", "pearson")
  ))


  pivot_assay_func <- switch(
    metric,
    "polarization" = PolarizationScoresToAssay,
    "colocalization" = ColocalizationScoresToAssay
  )
  spatial_metric_wide_format <- pivot_assay_func(spatial_metric_table, values_from, ...)

  # Make sure that the components match
  # Create an empty matrix (all 0's)
  tofillMat <- matrix(data = 0,
                      nrow = nrow(spatial_metric_wide_format),
                      ncol = ncol(object),
                      dimnames = list(rownames(spatial_metric_wide_format), colnames(object)))

  # Fill matrix where it overlaps
  # Any missing columns will be kept as 0's
  tofillMat[, colnames(spatial_metric_wide_format)] <- spatial_metric_wide_format

  return(tofillMat)
}
