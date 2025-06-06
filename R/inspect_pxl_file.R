#' Inspect a PXL file with MPX data
#'
#' This function inspects a PXL file and returns a tibble with
#' information about the files contained in the PXL file.
#'
#' @param pxl_file Path to a PXL file with MPX data
#'
#' @examples
#' pxl_file <- minimal_mpx_pxl_file()
#' inspect_pxl_file(pxl_file)
#'
#' @return A tibble with information about the files contained
#' in the PXL file
#'
#' @export
#'
inspect_pxl_file <- function(
  pxl_file
) {
  assert_file_exists(pxl_file)

  # List content of PXL file
  pxl_file_content <- utils::unzip(pxl_file, list = TRUE)

  adata_file <- grep("adata.h5ad", pxl_file_content$Name)
  edgelist_file <- grep("edgelist.parquet", pxl_file_content$Name)
  metadata_file <- grep("metadata.json", pxl_file_content$Name)
  polarization_file <- grep("polarization.parquet", pxl_file_content$Name)
  colocalization_file <- grep("colocalization.parquet", pxl_file_content$Name)
  layouts_dir <- grep("layouts.parquet", pxl_file_content$Name)

  if (length(layouts_dir) > 0) {
    start_marker <- "layout="
    end_marker <- "/component"
    layout_type <- stringr::str_extract(
      string = pxl_file_content$Name[layouts_dir],
      pattern = paste0(
        "(?<=\\b", start_marker,
        "\\b).*?(?=\\b", end_marker, "\\b)"
      )
    )
    start_marker <- "component="
    end_marker <- "/part"
    component <- stringr::str_extract(
      string = pxl_file_content$Name[layouts_dir],
      pattern = paste0("(?<=\\b", start_marker, "\\b).*?(?=\\b", end_marker, "\\b)")
    )
  }

  pxl_info <- tibble(
    file_type = c(
      "adata.h5ad",
      "edgelist.parquet",
      "metadata.json",
      "polarization.parquet",
      "colocalization.parquet"
    ),
    n = c(
      length(adata_file),
      length(edgelist_file),
      length(metadata_file),
      length(polarization_file),
      length(colocalization_file)
    ),
    file = list(
      pxl_file_content$Name[adata_file],
      pxl_file_content$Name[edgelist_file],
      pxl_file_content$Name[metadata_file],
      pxl_file_content$Name[polarization_file],
      pxl_file_content$Name[colocalization_file]
    )
  )

  if (length(layouts_dir) > 0) {
    pxl_info <- bind_rows(
      pxl_info,
      tibble(
        file_type = "layouts.parquet",
        n = length(layouts_dir),
        file = list(
          tibble(
            file = pxl_file_content$Name[layouts_dir],
            component = component,
            type = layout_type
          )
        )
      )
    )
  }

  return(pxl_info)
}
