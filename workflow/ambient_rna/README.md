# Ambient RNA removal

This module removes ambient RNA (cell-free RNA from the cell suspension that is captured in every droplet) from the counts of called cells.
Ambient RNA removal needs the raw, unfiltered droplet matrix of each library (all barcodes, including empty droplets), since the empty droplets show what the ambient RNA looks like.
The module runs one job per library (batch) and writes the corrected counts back into the input AnnData object as a new layer, so that the cells, their metadata and the raw counts stay as they are.

Available methods:

* `cellbender`: [CellBender](https://cellbender.readthedocs.io) `remove-background` ([Fleming et al. 2023](https://doi.org/10.1038/s41592-023-01943-7)), a deep generative model of the raw droplet matrix that estimates ambient RNA and barcode swapping per droplet.

## Environments

The following environments are needed for the module. Install only the ones you need.

- [`scanpy`](https://github.com/HCA-integration/scAtlasTb/blob/main/envs/scanpy.yaml)
- [`cellbender`](https://github.com/HCA-integration/scAtlasTb/blob/main/envs/cellbender.yaml) for `cellbender`

> Note: CellBender is practically GPU-only: it takes about half an hour per library on a GPU and many hours on a CPU.
> Set `use_gpu: true` in the config and create the `cellbender` environment on a machine with a GPU, so that conda installs the CUDA build of PyTorch (or set `CONDA_OVERRIDE_CUDA`, e.g. `CONDA_OVERRIDE_CUDA=12.6 conda env create -f envs/cellbender.yaml`).

## Configuration

```yaml
DATASETS:
  test:
    input:
      ambient_rna:
        cells: test/input/cells.h5ad
    ambient_rna:
      counts: layers/counts
      batch: batch
      raw_file: raw_file
      barcode: barcode
      gene_id_column: gene_ids
      methods:
        cellbender:
          expected_cells: 300
          total_droplets_included: 2000
          epochs: 50
```

* `counts`: Slot in the AnnData object that contains raw counts, e.g. `X`, `raw/X` or `layers/<layer_name>`.
  Cells that CellBender has not analyzed keep these counts in the output layer.
* `batch`: Column in `.obs` with the library of each cell. Each library is corrected separately, so this must be the unit that the raw droplet matrices belong to (e.g. one 10x channel). If not set, all cells must come from the same raw matrix.
* `raw_file`: Column in `.obs` with the path to the raw droplet matrix of each cell's library, e.g. the CellRanger `raw_feature_bc_matrix.h5`. All cells of a library must have the same path. CellBender reads CellRanger `.h5` files and MTX directories, AnnData `.h5ad` files and other formats (see the [CellBender documentation](https://cellbender.readthedocs.io/en/latest/reference/index.html)); the matrix must contain all barcodes, not only the called cells.
* `barcode`: Column in `.obs` with each cell's barcode as it appears in the raw matrix (e.g. `AAACCCAAGAAACACT-1`). If not set, the index of `.obs` is used.
* `gene_id_column`: Column in `.var` with the feature IDs of the raw matrix (e.g. Ensembl IDs). If not set, the index of `.var` is used. Genes are matched to the IDs of the raw matrix, or else to its feature names after making them unique like `var_names_make_unique`. Every gene of the input must be found.
* `methods`: Ambient RNA removal methods to run, either a list of method names or a mapping of method names to their parameters. Default: `cellbender`.
  * `cellbender`: parameters of `cellbender remove-background`, with underscores or dashes (e.g. `epochs`, `expected_cells`, `total_droplets_included`, `fpr`, `learning_rate`; see `cellbender remove-background --help`). Use `true` for flags without value. `input`, `output` and `cuda` are set by the pipeline (`--cuda` when `use_gpu: true` and a GPU is available), and only one `fpr` value is supported.

## Output

* `<out_dir>/ambient_rna/dataset~<dataset>/file_id~<file_id>.zarr` — the input AnnData with
  * `.layers['cellbender']`: corrected counts for every cell and gene of the input;
  * `.obs` columns:
    * `cellbender_cell_probability`: posterior probability that the droplet contains a cell (CellBender calls cells at > 0.5);
    * `cellbender_background_fraction`: fraction of the droplet's counts removed as background, as CellBender reports it (0 for droplets that it does not call as cells);
    * `cellbender_cell_size`, `cellbender_droplet_efficiency`: CellBender's estimates of the cell's true size and of the droplet's capture efficiency;
    * `cellbender_analyzed`: whether CellBender analyzed the barcode. Cells that CellBender did not analyze (because they are not among the `total_droplets_included` barcodes with the most counts, or not in the raw matrix) have missing values in the columns above and keep their raw counts;
    * `cellbender_removed_fraction`: fraction of the cell's counts in the input genes that was removed.
* `<image_dir>/ambient_rna/dataset~<dataset>/file_id~<file_id>/summary.tsv` — per library: number of cells, cells not analyzed, cells that CellBender also calls, median cell probability and the fraction of counts removed from the analyzed cells.
* `<out_dir>/ambient_rna/scatter/dataset~<dataset>/file_id~<file_id>/cellbender/<batch>/` — CellBender's own output per library: full (`cellbender.h5`) and cell-only (`cellbender_filtered.h5`) count matrices, a summary report (`cellbender.pdf`, `cellbender_report.html`), metrics (`cellbender_metrics.csv`), the log and the checkpoint (`ckpt.tar.gz`), from which a rerun resumes.

> Note: CellBender writes no counts for droplets that it does not call as cells (`cellbender_cell_probability` ≤ 0.5), so such cells are all-zero in `.layers['cellbender']` (`cellbender_removed_fraction` = 1).
> These are typically barcodes with few counts that the upstream cell calling kept.
> Remove them (e.g. keep `cellbender_cell_probability > 0.5`) before using the corrected counts, and check the per-library summary and CellBender's reports.

To use the corrected counts downstream, point the next module to the layer, e.g. `counts: layers/cellbender` for QC or `raw_counts: layers/cellbender` for integration.

## Test

The test data in `test/input` are synthetic (`test/generate_test_data.py`): two libraries of 300 cells, whose counts are 10% ambient RNA, and 7,700 empty droplets.
The test runs CellBender on CPU:

```commandline
bash workflow/ambient_rna/test/run_test.sh -c 4
```
