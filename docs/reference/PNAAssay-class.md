# The PNAAssay class

The PNAAssay object is an extended
[`Assay`](https://satijalab.github.io/seurat-object/reference/Assay-class.html)
for storing PNA single-cell data.

## Details

Compared to the
[`Assay`](https://satijalab.github.io/seurat-object/reference/Assay-class.html)
class, the PNAAssay class has three additional slots:

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
