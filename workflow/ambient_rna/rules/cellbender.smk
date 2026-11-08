def get_gpu_resource(resource_key, attempt):
    # stay on the GPU profile for retries: CellBender on CPU takes hours per batch
    return mcfg.get_resource(profile='gpu', resource_key=resource_key, attempt=attempt, attempt_to_cpu=attempt)


def parse_table_value(value):
    """Read a table cell as an int, a float, a boolean or else a string"""
    for cast in (int, float):
        try:
            return cast(value)
        except ValueError:
            pass
    return {'true': True, 'false': False}.get(value.lower(), value)


def read_batch_args(path, batch):
    """
    Arguments for one batch from a table with a column `batch` and one column per argument.
    Empty cells, and batches missing from the table, keep the method's arguments.
    """
    import pandas as pd

    table = pd.read_table(path, dtype=str, keep_default_na=False)
    if 'batch' not in table.columns:
        raise ValueError(f'batch_args table {path} needs a column "batch"')
    if not table['batch'].is_unique:
        raise ValueError(f'batch_args table {path} lists a batch more than once')
    rows = table[table['batch'] == batch]
    if rows.empty:
        return {}
    return {key: parse_table_value(value) for key, value in rows.iloc[0].drop('batch').items() if value != ''}


def get_method_args(wildcards, method):
    """A method's arguments for one batch: the method's own, updated from its `batch_args` table if set"""
    args = dict(get_methods(wildcards).get(method, {}))
    batch_args = args.pop('batch_args', None)
    if batch_args:
        args.update(read_batch_args(batch_args, wildcards.batch))
    return args


rule cellbender:
    """
    Run CellBender remove-background on the raw droplet matrix of one batch
    """
    input:
        batch=get_batch_file,
    output:
        h5=mcfg.out_dir / 'scatter' / params.wildcard_pattern / 'cellbender' / '{batch}' / 'cellbender.h5',
    params:
        args=lambda wildcards: get_method_args(wildcards, 'cellbender'),
        use_gpu=get_use_gpu(config),
    conda:
        get_env(config, 'cellbender')
    threads: 4
    resources:
        partition=lambda w, attempt: get_gpu_resource('partition', attempt),
        qos=lambda w, attempt: get_gpu_resource('qos', attempt),
        gpu=lambda w, attempt: get_gpu_resource('gpu', attempt),
        mem_mb=lambda w, attempt: get_gpu_resource('mem_mb', attempt),
    script:
        '../scripts/cellbender.py'
