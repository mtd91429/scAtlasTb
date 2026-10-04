"""
Write one file per batch with the path of the batch's raw droplet matrix
"""
import logging
logging.basicConfig(level=logging.INFO)
import re
from pathlib import Path
import pandas as pd

from utils.io import read_anndata

input_file = snakemake.input.anndata
output_dir = Path(snakemake.output.batches)
batch_key = snakemake.params.get('batch_key')
raw_file_key = snakemake.params.get('raw_file_key')

obs = read_anndata(input_file, obs='obs').obs
assert raw_file_key in obs.columns, f'Raw file column "{raw_file_key}" not found in obs columns.'

if batch_key is None or batch_key == 'None':
    logging.info('No batch key specified, all cells are processed as one batch.')
    groups = [('no_batch', obs)]
else:
    assert batch_key in obs.columns, f'Batch key "{batch_key}" not found in obs columns.'
    groups = obs.groupby(batch_key, observed=True, dropna=True)
    n_missing = obs[batch_key].isna().sum()
    if n_missing > 0:
        logging.warning(f'{n_missing} cells without a batch are not processed.')

output_dir.mkdir(exist_ok=True, parents=True)
file_names = set()
for batch, df in groups:
    raw_files = df[raw_file_key].dropna().astype(str).unique()
    assert len(raw_files) == 1, \
        f'Batch "{batch}" must have exactly one raw file in obs["{raw_file_key}"], found: {raw_files}'

    # batch names are used as file names
    file_name = re.sub(r'[^A-Za-z0-9_.-]', '_', str(batch))
    assert file_name not in file_names, f'Batches map to the same file name "{file_name}", rename the batches.'
    file_names.add(file_name)

    batch_file = output_dir / f'{file_name}.tsv'
    logging.info(f'Batch "{batch}": {df.shape[0]} cells, raw file {raw_files[0]} -> {batch_file}')
    pd.DataFrame(
        {'batch': [str(batch)], 'raw_file': [raw_files[0]], 'n_cells': [df.shape[0]]}
    ).to_csv(batch_file, sep='\t', index=False)
