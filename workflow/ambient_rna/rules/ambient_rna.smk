AMBIENT_RNA_METHODS = ['cellbender']


def get_methods(wildcards):
    """
    Methods configured for a dataset, as a mapping of method name to its parameters.
    `methods` can be a list of method names or a mapping of method name to parameters.
    """
    methods = mcfg.get_from_parameters(wildcards, 'methods', default=['cellbender'], check_query_keys=False)
    if isinstance(methods, str):
        methods = [methods]
    if not isinstance(methods, dict):
        methods = {method: {} for method in methods}
    methods = {method: (args if args else {}) for method, args in methods.items()}
    unknown = [method for method in methods if method not in AMBIENT_RNA_METHODS]
    if unknown:
        raise ValueError(f'Unknown ambient RNA method(s) {unknown}, available: {AMBIENT_RNA_METHODS}')
    return methods


checkpoint split_raw_files:
    input:
        anndata=lambda wildcards: mcfg.get_input_file(**wildcards)
    output:
        batches=directory(mcfg.out_dir / 'scatter' / params.wildcard_pattern / '_batches'),
    params:
        batch_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'batch_key'),
        raw_file_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'raw_file', check_null=True),
    conda:
        get_env(config, 'scanpy')
    resources:
        mem_mb=mcfg.get_resource(profile='cpu', resource_key='mem_mb'),
    script:
        '../scripts/split_raw_files.py'


def get_batch_dir(wildcards):
    return checkpoints.split_raw_files.get(**wildcards).output.batches


def get_batch_file(wildcards):
    return f'{get_batch_dir(wildcards)}/{wildcards.batch}.tsv'


def get_batches(wildcards):
    return glob_wildcards(f'{get_batch_dir(wildcards)}/{{batch}}.tsv').batch


def expand_per_batch(pattern, wildcards):
    """Expand a per-batch output pattern for all batches of a file"""
    return expand(pattern, batch=get_batches(wildcards), **wildcards)
