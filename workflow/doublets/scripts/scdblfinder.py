"""
scDblFinder (https://bioconductor.org/packages/scDblFinder) on each batch of a group of batches

The method runs in R (scdblfinder.R): this script reads the raw counts of the group, passes them
with the batch of each cell to R through an HDF5 file and writes the doublet score and call of
every cell. Each batch is processed on its own (one scDblFinder() call per batch, not the
`samples` argument, which selects features over all samples), with the random seed set before
each batch, so the result of a batch does not depend on the other batches of its group.

Parameters (``methods: {scdblfinder: {...}}``):
    seed: random seed, set before each batch (default: 0)
    any other: passed to scDblFinder() under its R name (e.g. dbr, dbr.sd, clusters, nfeatures)

Batches with fewer than 100 cells and cells without counts are not analysed: they get a score
of 0 and the call 'singlet', as in the other doublet methods of this module.
"""
import json
import logging
logging.basicConfig(level=logging.INFO)
import os
import subprocess
import tempfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from utils.io import read_anndata
from utils.misc import dask_compute

input_zarr = snakemake.input.zarr
output_tsv = snakemake.output.tsv
batches_txt = snakemake.input.batch
layer = snakemake.params.get('layer', 'X')
batch_key = snakemake.params.get('batch_key')
args = snakemake.params.get('args') or {}
threads = snakemake.threads
tmpdir = snakemake.resources.get('tmpdir')

MIN_CELLS = 100
R_SCRIPT = Path(__file__).resolve().with_suffix('.R')

args = dict(args)
seed = args.pop('seed', 0)
# set by the pipeline
reserved = {'sce', 'samples', 'returnType'}
assert not reserved & set(args), f'Parameters {sorted(reserved & set(args))} are set by the pipeline, remove them'

if batch_key == 'None':
    batch_key = None

logging.info(f'Read {input_zarr}...')
adata = read_anndata(
    input_zarr,
    X=layer,
    obs='obs',
    backed=True,
    dask=True,
)

with open(batches_txt, 'r') as f:
    batches = [line.strip().split('\t')[0] for line in f if line.strip()]

if batch_key is not None:
    adata = adata[adata.obs[batch_key].astype(str).isin(batches), :].copy()
    logging.info(f'Subset to {adata.n_obs} cells with specified batches.')
adata = dask_compute(adata, layers='X')

if batch_key is None:
    batch = pd.Series('no_batch', index=adata.obs_names)
else:
    batch = adata.obs[batch_key].astype(str)

X = sparse.csr_matrix(adata.X, dtype='float64')
X.sum_duplicates()
X.sort_indices()
sample = X.data[:100_000]
if not np.array_equal(sample, np.round(sample)):
    logging.warning(f'Counts in "{layer}" are not integers: scDblFinder expects raw counts')

# cells to analyse: cells with counts, in batches with at least MIN_CELLS of them
has_counts = np.asarray(X.sum(axis=1)).ravel() > 0
if not has_counts.all():
    logging.warning(f'{(~has_counts).sum()} cells without counts are not analysed')
batch_size = batch[has_counts].value_counts()
small_batches = batch_size.index[batch_size < MIN_CELLS]
if len(small_batches) > 0:
    logging.warning(f'Batches with fewer than {MIN_CELLS} cells with counts are not analysed: {small_batches.tolist()}')
analysed = has_counts & ~batch.isin(small_batches).to_numpy()

df = pd.DataFrame(
    {'scdblfinder_score': 0.0, 'scdblfinder_prediction': 'singlet'},
    index=adata.obs_names,
)

if analysed.any():
    X = X[analysed]
    if X.nnz > np.iinfo(np.int32).max:
        raise ValueError(
            f'{X.nnz} non-zero counts exceed the 2^31 - 1 entries that a sparse R matrix '
            '(dgCMatrix) can hold. Use a smaller chunk_size.'
        )
    batch_codes, batch_names = pd.factorize(batch[analysed], sort=True)
    var_names = adata.var_names.astype(str).to_numpy(dtype=object)
    del adata

    with tempfile.TemporaryDirectory(prefix='scdblfinder_', dir=tmpdir) as workdir:
        workdir = Path(workdir).resolve()
        input_file = workdir / 'input.h5'
        output_file = workdir / 'output.h5'

        # cells x genes CSR = genes x cells CSC as read by R
        logging.info(f'Write R input to {input_file}...')
        string_dtype = h5py.string_dtype()
        with h5py.File(input_file, 'w') as f:
            counts = f.create_group('counts')
            counts.create_dataset('data', data=X.data)
            counts.create_dataset('indices', data=X.indices.astype(np.int32))
            counts.create_dataset('indptr', data=X.indptr.astype(np.int32))
            counts.create_dataset('shape', data=np.array(X.shape, dtype=np.int32))
            f.create_dataset('var_names', data=var_names, dtype=string_dtype)
            f.create_dataset('batch', data=batch_codes.astype(np.int32))
            f.create_dataset('batch_names', data=batch_names.to_numpy(dtype=object), dtype=string_dtype)
        n_cells = X.shape[0]
        del X

        params = {
            'input': str(input_file),
            'output': str(output_file),
            'seed': seed,
            'args': args,
        }
        params_file = workdir / 'params.json'
        params_file.write_text(json.dumps(params, indent=2))

        env = os.environ.copy()
        for var in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']:
            env[var] = str(threads)
        command = ['Rscript', '--vanilla', str(R_SCRIPT), str(params_file)]
        logging.info(f'Run {" ".join(command)} with parameters:\n{json.dumps(params, indent=2)}')
        subprocess.run(command, check=True, env=env, cwd=workdir)

        with h5py.File(output_file, 'r') as f:
            score = f['score'][:]
            doublet = f['doublet'][:]
    assert score.shape[0] == n_cells and doublet.shape[0] == n_cells, \
        f'R returned {score.shape[0]} scores and {doublet.shape[0]} calls for {n_cells} cells'

    df.loc[analysed, 'scdblfinder_score'] = score
    df.loc[analysed, 'scdblfinder_prediction'] = np.where(doublet == 1, 'doublet', 'singlet')
else:
    logging.warning('No cells to analyse')

logging.info(f'{(df["scdblfinder_prediction"] == "doublet").sum()} of {df.shape[0]} cells called doublets')
df.to_csv(output_tsv, sep='\t')
