# Doublet detection

This module runs per-batch doublet detection using one or more callers (currently `scrublet`, `doubletdetection` and `scdblfinder`). The pipeline runs doublet callers within each library/batch and writes scores and predictions back into the AnnData object.

## Environments

The following environments are useful for running the module. Install only the ones you need.

- [`scanpy`](https://github.com/HCA-integration/scAtlasTb/blob/main/envs/scanpy.yaml)
- [`qc`](https://github.com/HCA-integration/scAtlasTb/blob/main/envs/qc.yaml) for `scrublet` (or `rapids_singlecell` with `use_gpu: true`) and `doubletdetection`
- [`scdblfinder`](https://github.com/HCA-integration/scAtlasTb/blob/main/envs/scdblfinder.yaml) for `scdblfinder` (Python, R 4.5 and [scDblFinder](https://bioconductor.org/packages/scDblFinder) 1.24 from Bioconductor)

# Configuration

```yaml
DATASETS:
  Lee2020:
    input:
      doublets: 
        Lee2020: test/input/load_data/harmonize_metadata/Lee2020.zarr
    doublets:
      counts: X
      batch: donor
      chunk_size: 10_000

  test:
    input:
      doublets:
        test: test/input/pbmc68k.h5ad
        test2: test/input/pbmc68k.h5ad
    doublets:
      counts: layers/counts
      methods:
        - scrublet
        - doubletdetection
        - scdblfinder

  test_scdblfinder:
    input:
      doublets:
        test: test/input/pbmc68k.h5ad
    doublets:
      counts: layers/counts
      methods:
        scrublet:
        scdblfinder:
          seed: 1
          dbr.per1k: 0.008  # scDblFinder() argument under its R name: expected doublet rate per 1,000 cells

defaults:
  datasets:
    - test
    - test_scdblfinder
    - Lee2020
```

* `counts`: Slot in anndata that contains raw (unnormalized) counts. Examples: `X`, `raw/X`, or `layers/<layer_name>`.
* `batch`: Column in `obs` that contains batch/library IDs. Doublet detection is executed separately per batch.
* `chunk_size`: Number of cells used to group batches for more efficient parallel processing. Default: `100_000`.
* `methods`: Doublet callers to run for this dataset, either a list of method names (e.g. `['scrublet', 'doubletdetection']`) or a mapping of method names to their parameters. Available: `scrublet`, `doubletdetection`, `scdblfinder`. Default: `scrublet`.
  Parameters are only supported for `scdblfinder`:
  * `seed`: random seed, set before each batch. Default: `0`.
  * any other key is passed to [`scDblFinder()`](https://bioconductor.org/packages/scDblFinder) under its R name (e.g. `dbr`, `dbr.sd`, `clusters`, `nfeatures`), with scDblFinder's defaults otherwise (random artificial doublets, expected doublet rate of 0.8% per 1,000 cells).

> Note: `counts` is resolved relative to the input object (e.g., anndata.X or anndata.layers). Methods are run per-batch; choose `batch` to reflect library-level grouping so cross-library doublets are not considered.

> `scdblfinder` runs [scDblFinder](https://bioconductor.org/packages/scDblFinder) in R, one `scDblFinder()` call per batch with the random seed set before each call, so the result of a batch does not depend on how batches are grouped into chunks. It does not use scDblFinder's `samples` argument, which would select the features over all batches of a chunk. Batches with fewer than 100 cells, and cells without counts, are not analysed: as with the other callers, they get a score of 0 and are called singlets.

## Output

Results are written into the AnnData `.obs` and to summary files/plots:

* `<out_dir>/doublets/dataset~<datasets>/file_id~<file_id>.zarr` — Updated AnnData with caller-specific score and prediction columns in `.obs`, for example:
  - scrublet:
    - `scrublet_score`
    - `scrublet_prediction`
  - doubletdetection:
    - `doubletdetection_score`
    - `doubletdetection_prediction`
  - scdblfinder:
    - `scdblfinder_score`: scDblFinder's doublet score (`scDblFinder.score`)
    - `scdblfinder_prediction`: `doublet` or `singlet` (`scDblFinder.class`)
* `<image_dir>/doublets/dataset~<datasets>/file_id~<file_id>/doublet_summary.tsv` — Summary table with per-batch and per-method statistics.
* `<image_dir>/doublets/dataset~<datasets>/file_id~<file_id>/plots/` — Visualizations of score distributions and per-batch results (e.g., histograms, per-batch summary plots).

If additional callers are enabled, corresponding score/prediction columns will be added to `.obs` using method-specific names to avoid collisions.