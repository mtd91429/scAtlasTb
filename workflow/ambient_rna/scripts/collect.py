"""
Collect the ambient RNA-corrected counts of all batches into the input AnnData
"""
import logging
logging.basicConfig(level=logging.INFO)
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from anndata.utils import make_index_unique

from utils.io import read_anndata, write_zarr_linked


input_file = snakemake.input.anndata
batch_files = snakemake.input.batches
cellbender_files = snakemake.input.get('cellbender', [])
output_zarr = snakemake.output.zarr
output_summary = snakemake.output.summary
counts_slot = snakemake.params.get('counts', 'X')
batch_key = snakemake.params.get('batch_key')
barcode_key = snakemake.params.get('barcode_key')
gene_id_column = snakemake.params.get('gene_id_column')

LAYER = 'cellbender'
# per-droplet latent variables of CellBender, saved as obs columns
LATENTS = ['cell_probability', 'background_fraction', 'cell_size', 'droplet_efficiency']
CELL_PROBABILITY_CUTOFF = 0.5  # CellBender's cell call


def read_cellbender_h5(file):
    """Read corrected counts (barcodes x features) and latent variables from a CellBender output file"""
    with h5py.File(file, 'r') as f:
        group = f['matrix']
        n_features, n_barcodes = group['shape'][()]
        # stored as features x barcodes, like CellRanger
        matrix = sparse.csc_matrix(
            (group['data'][()], group['indices'][()], group['indptr'][()]),
            shape=(n_features, n_barcodes),
        ).T.tocsr()
        barcodes = group['barcodes'][()].astype(str)
        feature_ids = group['features/id'][()].astype(str)
        feature_names = group['features/name'][()].astype(str)
        latent_group = f['droplet_latents']
        latents = pd.DataFrame(
            {key: latent_group[key][()] for key in LATENTS if key in latent_group},
            index=latent_group['barcode_indices_for_latents'][()],
        ).reindex(columns=LATENTS)
    return matrix, barcodes, feature_ids, feature_names, latents


def match_genes(genes, feature_ids, feature_names):
    """Positions of the input genes among CellBender's features, matched by feature ID or by name"""
    candidates = {
        'feature IDs': pd.Index(feature_ids),
        # same names as scanpy.read_10x_h5 followed by var_names_make_unique
        'feature names': make_index_unique(pd.Index(feature_names)),
    }
    for label, features in candidates.items():
        if not features.is_unique:
            continue
        positions = features.get_indexer(genes)
        if (positions >= 0).all():
            logging.info(f'Matched {len(genes)} genes to CellBender {label}')
            return positions
    missing = genes[pd.Index(feature_ids).get_indexer(genes) < 0]
    raise ValueError(
        f'{len(missing)} of {len(genes)} genes are neither CellBender feature IDs nor names, e.g. {list(missing[:5])}. '
        'Set `gene_id_column` to a var column with the feature IDs of the raw matrix.'
    )


logging.info(f'Read {input_file}...')
if input_file.endswith('.h5ad'):
    adata = read_anndata(input_file)
else:
    adata = read_anndata(input_file, obs='obs', var='var')
counts = read_anndata(input_file, X=counts_slot).X
counts = sparse.csr_matrix(counts)
logging.info(adata.__str__())

genes = pd.Index(adata.var[gene_id_column] if gene_id_column else adata.var_names).astype(str)
barcodes = pd.Index(adata.obs[barcode_key] if barcode_key else adata.obs_names).astype(str)
if batch_key is None or batch_key == 'None':
    batches = pd.Series('no_batch', index=adata.obs_names)
else:
    batches = adata.obs[batch_key].astype(str)

cellbender_files = {Path(file).parent.name: file for file in cellbender_files}
latent_df = pd.DataFrame(np.nan, index=adata.obs_names, columns=LATENTS)
analyzed = np.zeros(adata.n_obs, dtype=bool)
corrected_blocks = []
summary = []

for batch_file in batch_files:
    batch_info = pd.read_table(batch_file, dtype=str).iloc[0]
    batch = batch_info['batch']
    cells = np.flatnonzero(batches.values == batch)
    file = cellbender_files[Path(batch_file).stem]
    logging.info(f'Batch "{batch}": {len(cells)} cells, read {file}...')
    matrix, cb_barcodes, feature_ids, feature_names, latents = read_cellbender_h5(file)

    # input cells among the droplets that CellBender analyzed
    positions = pd.Index(cb_barcodes).get_indexer(barcodes[cells])
    found = positions >= 0
    if not found.any():
        raise ValueError(
            f'None of the {len(cells)} barcodes of batch "{batch}" are in the raw matrix, '
            f'e.g. {list(barcodes[cells][:3])} vs. {list(cb_barcodes[:3])}. '
            'Set `barcode` to an obs column with barcodes as they appear in the raw matrix.'
        )
    is_analyzed = found & np.isin(positions, latents.index)
    cells_analyzed = cells[is_analyzed]
    analyzed[cells_analyzed] = True
    latent_df.iloc[cells_analyzed] = latents.loc[positions[is_analyzed], LATENTS].values

    gene_positions = match_genes(genes, feature_ids, feature_names)
    block = matrix[positions[is_analyzed]][:, gene_positions].astype(counts.dtype)
    corrected_blocks.append((cells_analyzed, block))

    raw_total = counts[cells_analyzed].sum()
    cell_probability = latent_df['cell_probability'].values[cells_analyzed]
    summary.append({
        'batch': batch,
        'raw_file': batch_info['raw_file'],
        'n_cells': len(cells),
        'n_cells_not_in_raw_file': int((~found).sum()),
        'n_cells_not_analyzed': int((found & ~is_analyzed).sum()),
        'n_cells_called': int((cell_probability > CELL_PROBABILITY_CUTOFF).sum()),
        'median_cell_probability': np.median(cell_probability) if len(cell_probability) else np.nan,
        'counts_raw': raw_total,
        'counts_corrected': block.sum(),
        'fraction_removed': 1 - block.sum() / raw_total if raw_total > 0 else np.nan,
    })

# cells that CellBender did not analyze keep their raw counts
not_analyzed = np.flatnonzero(~analyzed)
logging.info(f'{analyzed.sum()} cells corrected, {len(not_analyzed)} cells keep their raw counts')
order = np.concatenate([cells for cells, _ in corrected_blocks] + [not_analyzed])
corrected = sparse.vstack(
    [block for _, block in corrected_blocks] + [counts[not_analyzed]],
    format='csr',
)[np.argsort(order)]
assert corrected.shape == counts.shape

raw_total = np.asarray(counts.sum(axis=1)).ravel()
corrected_total = np.asarray(corrected.sum(axis=1)).ravel()
with np.errstate(divide='ignore', invalid='ignore'):
    removed_fraction = np.where(raw_total > 0, 1 - corrected_total / raw_total, np.nan)

for column in LATENTS:
    adata.obs[f'cellbender_{column}'] = latent_df[column].values
adata.obs['cellbender_analyzed'] = analyzed
adata.obs['cellbender_removed_fraction'] = removed_fraction
adata.layers[LAYER] = corrected

summary = pd.DataFrame(summary, columns=[
    'batch', 'raw_file', 'n_cells', 'n_cells_not_in_raw_file', 'n_cells_not_analyzed',
    'n_cells_called', 'median_cell_probability', 'counts_raw', 'counts_corrected', 'fraction_removed',
])
logging.info(f'Summary per batch:\n{summary.to_string()}')
Path(output_summary).parent.mkdir(parents=True, exist_ok=True)
summary.to_csv(output_summary, sep='\t', index=False)

logging.info(f'Write {output_zarr}...')
logging.info(adata.__str__())
if input_file.endswith('.h5ad'):
    input_layers = []
else:
    layers_dir = Path(input_file) / 'layers'
    input_layers = [f.name for f in layers_dir.iterdir() if f.is_dir()] if layers_dir.is_dir() else []
write_zarr_linked(
    adata,
    input_file,
    output_zarr,
    files_to_keep=['obs', 'var', f'layers/{LAYER}'],
    slot_map={f'layers/{layer}': f'layers/{layer}' for layer in input_layers if layer != LAYER},
)
