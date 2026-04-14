# Export Seurat object MPX data to a PXL file

**\[experimental\]**

The function only exports the essential data required to create a
functional PXL file. See details below for a description of what data is
exported.

## Usage

``` r
WriteMPX_pxl_file(
  object,
  file,
  assay = NULL,
  export_layouts = FALSE,
  overwrite = FALSE
)
```

## Arguments

- object:

  A `Seurat` object with a `CellGraphAssay5` assay object created with
  pixelatorR.

- file:

  A character string specifying the path to the PXL file to be created.

- assay:

  A character string specifying the name of the `CellGraphAssay5`. If
  set to `NULL` the default assay will be used.

- export_layouts:

  A logical value specifying whether to export the layouts from the
  `Seurat` object. If set to `TRUE`, each component must have a layout.
  See [`LoadCellGraphs`](LoadCellGraphs.md) and
  [`ComputeLayout`](ComputeLayout.md) for details about how to load
  graphs and compute layouts.

- overwrite:

  A logical value specifying whether to overwrite the `file` if it
  already exists.

## Value

Nothing. The function writes the PXL file to the specified location.

## Exported data

- Count data and metadata : The raw count matrix and component metadata
  from the Seurat object are exported into a .h5ad file that can be read
  with the anndata Python library.

- Polarization scores : The polarization scores are exported to a
  .parquet file.

- Colocalization scores : The colocalization scores are exported to a
  .parquet file.

- Edgelist data : The edgelist data is first collected from the original
  PXL file(s), filtered to include the components currently available in
  the `Seurat` object, and then the component IDs are updated. The
  resulting merged edgelist data is exported to a .parquet file.

- Sample meta data : The sample meta data is extracted from the original
  PXL file(s), then merged and exported to a .json file.

The structure of the PXL file is detailed below:

\|– adata.h5ad  
\|– polarization.parquet  
\|– colocalization.parquet  
\|– metadata.json  
\|– edgelist.parquet

The merged files are converted into a zip archive and saved to the
target PXL file.

NOTE: Factors are currently not supported. These will be converted to
string arrays.

## Examples

``` r
# Use Assay5 as the default assay version
options(Seurat.object.assay.version = "v5")

# Create Seurat object
pxl_file <- minimal_mpx_pxl_file()
se <- ReadMPX_Seurat(pxl_file)
#> ! Failed to remove temporary file C:/Users/max/AppData/Local/Temp/RtmpmOhqBt/file62e82d5e2066.h5ad

se_merged <- merge(se, list(se, se, se))
#> Warning: Some cell names are duplicated across objects provided. Renaming to enforce unique cell names.
pxl_file <- fs::path_temp("small.pxl")
if (fs::file_exists(pxl_file)) {
  fs::file_delete(pxl_file)
}

# Export data to a new PXL file
WriteMPX_pxl_file(se_merged, pxl_file)
#> ── Creating PXL file ───────────────────────────────────────────────────────────
#> → Writing spatial metrics to .parquet files...
#> → Writing polarization scores to polarization.parquet
#> → Writing colocalization scores to colocalization.parquet
#> ✔ Exported spatial metrics
#> → Collecting edge lists from 4 files
#> → Extracting sample 1 edge list
#> → Extracting sample 2 edge list
#> → Extracting sample 3 edge list
#> → Extracting sample 4 edge list
#> →  Merging edge lists
#> → Loaded sample 1 edge list
#> → Loaded sample 2 edge list
#> → Loaded sample 3 edge list
#> → Loaded sample 4 edge list
#> → Merging edge lists
#> → Exporting merged edge list
#> ✔ Exported merged edge list
#> ✔ Exported merged meta data
#> ✔ Exported anndata file
#> ℹ Saving PXL file to C:/Users/max/AppData/Local/Temp/RtmpmOhqBt/small.pxl
#> ✔ Finished!

# Read the new PXL file
se_merged <- ReadMPX_Seurat(pxl_file)
#> ! Failed to remove temporary file C:/Users/max/AppData/Local/Temp/RtmpmOhqBt/file62e87eeb34a4.h5ad

# Reset global option
options(Seurat.object.assay.version = "v3")
```
