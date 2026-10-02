"""
STACAS integration (https://github.com/carmonalab/STACAS): Seurat integration anchors between
pairs of batches, rescored by STACAS and, optionally, filtered by cell type labels
(semi-supervised), then integrated along a guide tree of the batches.

The method runs in R (stacas.R, called through seurat_utils). STACAS is not available from
conda channels: envs/stacas.post-deploy.sh installs it into the stacas environment.
The method starts from the raw counts (layers/counts): each cell is log-normalized on all genes
(Seurat's LogNormalize) before the genes are subset. Output: embedding in obsm['X_emb'], the PCA
that STACAS computes on the integrated data.

Hyperparameters:
    anchor_features: genes used to find anchors. Default (none): the integration features
        (var_mask). A number lets STACAS select that many genes itself from all genes, after
        removing its gene block list (genesBlockList, default "default": cell cycle,
        mitochondrial, ribosomal and other genes). STACAS (2.4.1) applies the block list only
        in this mode.
    cell_labels: obs column with cell type labels for semi-supervised integration; cells with
        missing labels count as unlabelled. Default (none): unsupervised.
    dims: number of dimensions used to find anchors and of the returned embedding (default: 30)
    reference: batch, or list of batches, used as reference (default: none); pass a list inside
        a list, since lists are expanded into separate runs, e.g. reference: [[batch1, batch2]]
    scale_factor: LogNormalize scale factor (default: 10000)
    seed: random seed of STACAS, used only with label_confidence < 1 (default: the module seed)
    any other: passed to Run.STACAS() under its R name, with underscores replaced by dots
        (e.g. k_anchor: 5 -> k.anchor = 5, min_sample_size: 50 -> min.sample.size = 50)
"""
import logging
logging.basicConfig(level=logging.INFO)
from pathlib import Path
from pprint import pformat

from integration_utils import add_metadata, remove_slots, clean_categorical_column
from seurat_utils import read_counts, run_r_method, r_arguments
from utils.io import write_zarr_linked


input_file = snakemake.input[0]
output_file = snakemake.output[0]
wildcards = snakemake.wildcards
params = snakemake.params
batch_key = wildcards.batch

hyperparams = params.get('hyperparams', {})
hyperparams = {} if hyperparams is None else dict(hyperparams)

anchor_features = hyperparams.pop('anchor_features', None)
cell_labels = hyperparams.pop('cell_labels', None)
reference = hyperparams.pop('reference', None)
if reference is not None and not isinstance(reference, list):
    reference = [reference]
r_params = {
    'seed': hyperparams.pop('seed', params.get('seed', 0)),
    'dims': hyperparams.pop('dims', 30),
    'scale_factor': hyperparams.pop('scale_factor', 1e4),
    'anchor_features': None if anchor_features is None else int(anchor_features),
    'reference': None if reference is None else [str(batch) for batch in reference],
    'semi_supervised': cell_labels is not None,
    'stacas_args': r_arguments(hyperparams),
}

# STACAS selecting its own anchor features needs every gene
adata, library_size = read_counts(
    input_file,
    var_column='integration_features' if anchor_features is None else None,
)
clean_categorical_column(adata, batch_key)
obs = {'batch': adata.obs[batch_key]}
if cell_labels is not None:
    assert cell_labels in adata.obs.columns, f'Column {cell_labels} is missing'
    obs['cell_labels'] = adata.obs[cell_labels]

# run method
logging.info(f'Run STACAS with parameters {pformat(r_params)}...')
adata.obsm['X_emb'] = run_r_method(
    script=Path(__file__).with_suffix('.R'),
    adata=adata,
    library_size=library_size,
    obs=obs,
    params=r_params,
    threads=snakemake.threads,
    tmpdir=snakemake.resources.get('tmpdir'),
)

# prepare output adata
adata = remove_slots(adata=adata, output_type=params['output_type'])
add_metadata(adata, wildcards, params)

logging.info(f'Write {output_file}...')
logging.info(adata.__str__())
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=['obsm', 'uns'],
)
