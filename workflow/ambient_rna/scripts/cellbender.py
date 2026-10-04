"""
Run CellBender remove-background on the raw droplet matrix of one batch
"""
import logging
logging.basicConfig(level=logging.INFO)
import shlex
import subprocess
from pathlib import Path
import pandas as pd
import torch

batch_file = snakemake.input.batch
output_h5 = Path(snakemake.output.h5).resolve()
args = snakemake.params.get('args') or {}
use_gpu = snakemake.params.get('use_gpu', False)
threads = snakemake.threads

batch_info = pd.read_table(batch_file).iloc[0]
raw_file = Path(batch_info['raw_file']).resolve()
assert raw_file.exists(), f'Raw file for batch "{batch_info["batch"]}" not found: {raw_file}'

# set by the pipeline
reserved = {'input', 'output', 'cuda'}
args = {key.replace('-', '_'): value for key, value in args.items()}
assert not reserved & set(args), f'Parameters {sorted(reserved & set(args))} are set by the pipeline, remove them'
fpr = args.get('fpr')
if isinstance(fpr, (list, tuple)):
    assert len(fpr) == 1, f'Only one fpr value is supported, got {fpr}'
args.setdefault('cpu_threads', threads)

use_cuda = use_gpu and torch.cuda.is_available()
if use_gpu and not use_cuda:
    logging.warning('use_gpu is set, but no GPU is available: running CellBender on CPU (slow)')

cmd = [
    'cellbender', 'remove-background',
    '--input', str(raw_file),
    '--output', str(output_h5),
]
if use_cuda:
    cmd.append('--cuda')
for key, value in args.items():
    flag = f'--{key.replace("_", "-")}'
    if value is None or value is False:
        continue
    if value is True:
        cmd.append(flag)
    elif isinstance(value, (list, tuple)):
        cmd.extend([flag, *map(str, value)])
    else:
        cmd.extend([flag, str(value)])

# CellBender writes its checkpoint (ckpt.tar.gz) to the working directory, from which a rerun resumes
output_h5.parent.mkdir(parents=True, exist_ok=True)
logging.info(f'Batch "{batch_info["batch"]}" ({batch_info["n_cells"]} cells in the input)')
logging.info(f'Run in {output_h5.parent}:\n{shlex.join(cmd)}')
subprocess.run(cmd, check=True, cwd=output_h5.parent)
