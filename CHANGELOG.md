# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [UNRELEASED] - 2024-??-??

## [0.2.0] - 2024-03-18

### Added new test data set - five_cells.pxl

Once `pixelatorR` is installed, the test data can be read from the following path:

````
pxl_file <- system.file("extdata/five_cells",
  "five_cells.pxl",
  package = "pixelatorR"
)
````

- The data set contains 5 cells from three different MPX data sets:
   - [CD3 capped cell: RCVCMP0000217](https://software.pixelgen.com/datasets/) from the CD3 capping data set (stimulated)
   - [B cell: RCVCMP0000118, CD4 T cell: RCVCMP0000487, CD8 T cell: RCVCMP0000655](https://software.pixelgen.com/datasets/) from a PBMC data set
   - [Uropod cell: RCVCMP0000263](https://software.pixelgen.com/datasets/) from the RANTES stimulated data set

- All cells now have fully connected components
- The total file size is 2.1 MB, where the following fields have been excluded from the edgelist.parquet file: "umi", "sequence", "umi_unique_count", "upi_unique_count"

## [0.1.3] - 2024-03-18

- The latest version of `arrow` on CRAN (v15.0.1) contains several bugs which are picked up by the R CMD checks. Fixed the version of `arrow` to v14.0.2.1 until a stable release is available. 

## [0.1.2] - 2024-03-14

* Fixed bug in `LoadCellGraphs` when working on Seurat objects with multiple data sets (merged). Previous version returned components in an invalid order.

## [0.1.1] - 2024-02-29

* Fixed bug in `ComputeLayout` when normalizing or projecting layout coordinates
* Added layout utility functions:
  * `center_layout_coordinates` - Center layout coordinates at the origin
  * `normalize_layout_coordinates` - Normalize layout coordinates such that they are centered at origin and the have a median length (euclidean norm) of 1
  * `project_layout_coordinates_on_unit_sphere` - Project layout 3D coordinates onto the unit sphere

## [0.1.0] - 2024-02-21

* First public release of pixelatorR.
