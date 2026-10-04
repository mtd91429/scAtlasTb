"""
Generate synthetic test data for the ambient_rna module

Per batch, a raw droplet matrix in CellRanger v3 HDF5 format (all barcodes): 300 cells of 3 cell
types whose counts are 10% ambient RNA, and 7,700 empty droplets that contain only ambient RNA.
The ambient profile is the mixture of the cell type profiles. The cells are also saved as an
AnnData object (the "called cells" that the module corrects), with the batch, barcode and path to
the raw matrix in .obs.

Run from the module directory (workflow/ambient_rna) in the scanpy environment:
    python test/generate_test_data.py
"""
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

OUT_DIR = Path('test/input')
N_GENES = 200
N_CELL_TYPES = 3
N_CELLS = 300
N_EMPTY = 7_700
AMBIENT_FRACTION = 0.1


def write_10x_h5(file, matrix, barcodes, gene_ids, gene_names):
    """Write a barcodes x genes count matrix in CellRanger v3 HDF5 format"""
    matrix = sparse.csc_matrix(matrix.T)  # stored as genes x barcodes
    with h5py.File(file, 'w') as f:
        group = f.create_group('matrix')
        group.create_dataset('barcodes', data=np.array(barcodes, dtype='S'))
        group.create_dataset('data', data=matrix.data.astype(np.int32), compression='gzip')
        group.create_dataset('indices', data=matrix.indices.astype(np.int64), compression='gzip')
        group.create_dataset('indptr', data=matrix.indptr.astype(np.int64))
        group.create_dataset('shape', data=np.array(matrix.shape, dtype=np.int32))
        features = group.create_group('features')
        features.create_dataset('id', data=np.array(gene_ids, dtype='S'))
        features.create_dataset('name', data=np.array(gene_names, dtype='S'))
        features.create_dataset('feature_type', data=np.array(['Gene Expression'] * len(gene_ids), dtype='S'))
        features.create_dataset('genome', data=np.array(['test'] * len(gene_ids), dtype='S'))


def random_barcodes(rng, n):
    barcodes = set()
    while len(barcodes) < n:
        barcodes.add(''.join(rng.choice(list('ACGT'), 16)) + '-1')
    return sorted(barcodes)


rng = np.random.default_rng(0)
gene_ids = [f'TESTG{i:05d}' for i in range(N_GENES)]
gene_names = [f'Gene{i}' for i in range(N_GENES)]
profiles = rng.dirichlet(np.full(N_GENES, 0.1), size=N_CELL_TYPES)

OUT_DIR.mkdir(parents=True, exist_ok=True)
cells = []
for batch_number in (1, 2):
    batch = f'batch{batch_number}'
    cell_types = rng.choice(N_CELL_TYPES, size=N_CELLS, p=[0.5, 0.3, 0.2])
    ambient = np.bincount(cell_types, minlength=N_CELL_TYPES) @ profiles / N_CELLS

    cell_size = rng.lognormal(np.log(2_000), 0.4, size=N_CELLS)
    cell_counts = (
        rng.poisson(np.outer(cell_size * (1 - AMBIENT_FRACTION), np.ones(N_GENES)) * profiles[cell_types])
        + rng.poisson(np.outer(cell_size * AMBIENT_FRACTION, ambient))
    )
    empty_size = rng.lognormal(np.log(15), 0.5, size=N_EMPTY)
    empty_counts = rng.poisson(np.outer(empty_size, ambient))

    barcodes = random_barcodes(rng, N_CELLS + N_EMPTY)
    order = rng.permutation(N_CELLS + N_EMPTY)  # cells are not the first barcodes
    raw = sparse.csr_matrix(np.vstack([cell_counts, empty_counts])[np.argsort(order)])
    barcodes = np.array(barcodes)
    cell_barcodes = barcodes[order[:N_CELLS]]

    raw_file = OUT_DIR / f'raw_feature_bc_matrix_{batch}.h5'
    write_10x_h5(raw_file, raw, barcodes, gene_ids, gene_names)

    adata = ad.AnnData(
        X=sparse.csr_matrix(cell_counts, dtype=np.float32),
        obs=pd.DataFrame(
            {
                'batch': batch,
                'barcode': cell_barcodes,
                'raw_file': str(raw_file),
                'cell_type': [f'type{t}' for t in cell_types],
            },
            index=[f'{batch}:{bc}' for bc in cell_barcodes],
        ),
        var=pd.DataFrame({'gene_ids': gene_ids}, index=gene_names),
    )
    cells.append(adata)
    print(f'{batch}: {N_CELLS} cells, {N_EMPTY} empty droplets -> {raw_file}')

adata = ad.concat(cells, merge='same')
adata.layers['counts'] = adata.X.copy()
for column in ['batch', 'raw_file', 'cell_type']:
    adata.obs[column] = adata.obs[column].astype('category')
adata.write_h5ad(OUT_DIR / 'cells.h5ad', compression='gzip')
print(adata)
