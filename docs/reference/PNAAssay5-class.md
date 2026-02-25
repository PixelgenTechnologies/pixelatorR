# The PNAAssay5 class

The PNAAssay5 object is an extended
[`Assay5`](https://satijalab.github.io/seurat-object/reference/Assay5-class.html)
for storing PNA single-cell data.

## Details

Compared to the
[`Assay5`](https://satijalab.github.io/seurat-object/reference/Assay5-class.html)
class, the PNAAssay5 class has three additional slots:

- `cellgraphs`: A named list of [`CellGraph`](CellGraph-class.md)
  objects

- `proximity`: A `tbl_df` with proximity scores

- `fs_map`: A `tbl_df` with information on source PXL file

## Slots

- `cellgraphs`:

  A named list of [`CellGraph`](CellGraph-class.md) objects

- `proximity`:

  A `tbl_df` with proximity scores

- `fs_map`:

  A `tbl_df` with information on source PXL file paths, sample IDs, and
  component IDs
