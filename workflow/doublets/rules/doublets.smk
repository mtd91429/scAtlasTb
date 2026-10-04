DOUBLET_METHODS = ['scrublet', 'doubletdetection', 'scdblfinder']


def get_methods(wildcards):
    """
    Doublet methods configured for a dataset, as a mapping of method name to its parameters.
    `methods` can be a list of method names or a mapping of method names to parameters
    (parameters are only supported for scdblfinder).
    """
    methods = mcfg.get_from_parameters(wildcards, 'methods', default=['scrublet'])
    if isinstance(methods, str):
        methods = [methods]
    if not isinstance(methods, dict):
        methods = {method: {} for method in methods}
    methods = {method: (args if args else {}) for method, args in methods.items()}
    unknown = [method for method in methods if method not in DOUBLET_METHODS]
    if unknown:
        raise ValueError(f'Unknown doublet method(s) {unknown}, available: {DOUBLET_METHODS}')
    with_args = [method for method, args in methods.items() if args and method != 'scdblfinder']
    if with_args:
        raise ValueError(f'Parameters are only supported for scdblfinder, not for {with_args}')
    return methods


checkpoint split_batches:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards)
    output:
        batches=directory(mcfg.out_dir / 'scatter' /  params.wildcard_pattern / '_batches'),
    params:
        batch_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'batch_key'),
        chunk_size=lambda wildcards: mcfg.get_from_parameters(wildcards, 'chunk_size', default=100_000),
    conda:
        get_env(config, 'scanpy')
    script:
        '../scripts/split_batches.py'


def get_checkpoint_output(wildcards):
    return f'{checkpoints.split_batches.get(**wildcards).output[0]}/{{batch}}.txt'


def get_from_checkpoint(wildcards, pattern=None):
    checkpoint_output = get_checkpoint_output(wildcards)
    if pattern is None:
        pattern = checkpoint_output
    return expand(
        pattern,
        batch=glob_wildcards(checkpoint_output).batch,
        allow_missing=True
    )


def get_mem_mb(attempt, profile, factor=2):
    mem_mb = mcfg.get_resource(profile=profile, resource_key='mem_mb', attempt=attempt)
    try:
        mem_mb = int(mem_mb)
    except ValueError:
        return mem_mb
    return max(1_000, int(mem_mb // factor))


rule scrublet:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
        batch=get_checkpoint_output,
    output:
        tsv=mcfg.out_dir / 'scatter' / params.wildcard_pattern / 'scrublet' / '{batch}.tsv',
    params:
        batch_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'batch_key', check_query_keys=False),
        layer=lambda wildcards: mcfg.get_from_parameters(wildcards, 'counts', default='X'),
    conda:
        lambda wildcards: get_env(
            config,
            'qc',
            gpu_env='rapids_singlecell',
            no_gpu=not mcfg.get_from_parameters(wildcards, 'use_gpu', default=False)
        )
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='partition', attempt=attempt),
        qos=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='qos', attempt=attempt),
        gpu=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='gpu', attempt=attempt),
        mem_mb=lambda w, attempt: get_mem_mb(attempt=attempt, profile='gpu'),
    script:
        '../scripts/scrublet.py'


rule doubletdetection:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
        batch=get_checkpoint_output
    output:
        tsv=mcfg.out_dir / 'scatter' / params.wildcard_pattern / 'doubletdetection' / '{batch}.tsv',
    params:
        batch_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'batch_key', check_query_keys=False),
        layer=lambda wildcards: mcfg.get_from_parameters(wildcards, 'counts', default='X'),
    conda:
        get_env(config, 'qc')
    threads: 3
    resources:
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),
        gpu=mcfg.get_resource(profile='cpu',resource_key='gpu'),
        mem_mb=lambda w, attempt: get_mem_mb(attempt, profile='cpu'),
    # shadow: "minimal"
    script:
        '../scripts/doubletdetection.py'


rule scdblfinder:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
        batch=get_checkpoint_output
    output:
        tsv=mcfg.out_dir / 'scatter' / params.wildcard_pattern / 'scdblfinder' / '{batch}.tsv',
    params:
        batch_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'batch_key', check_query_keys=False),
        layer=lambda wildcards: mcfg.get_from_parameters(wildcards, 'counts', default='X'),
        args=lambda wildcards: get_methods(wildcards).get('scdblfinder', {}),
    conda:
        get_env(config, 'scdblfinder')
    resources:
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),
        gpu=mcfg.get_resource(profile='cpu',resource_key='gpu'),
        mem_mb=lambda w, attempt: get_mem_mb(attempt, profile='cpu'),
    script:
        '../scripts/scdblfinder.py'


def collect_results(wildcards):
    methods = get_methods(wildcards)
    files = {'zarr': mcfg.get_input_file(**wildcards)}
    if 'scrublet' in methods:
        files['scrublet'] = get_from_checkpoint(wildcards, rules.scrublet.output.tsv)
    if 'doubletdetection' in methods:
        files['doubletdetection'] = get_from_checkpoint(wildcards, rules.doubletdetection.output.tsv)
    if 'scdblfinder' in methods:
        files['scdblfinder'] = get_from_checkpoint(wildcards, rules.scdblfinder.output.tsv)
    return files


rule collect:
    input:
        unpack(collect_results)
        # zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
        # scrublet=lambda wildcards: get_from_checkpoint(wildcards, rules.scrublet.output.tsv),
        # doubletdetection=lambda wildcards: get_from_checkpoint(wildcards, rules.doubletdetection.output.tsv),
    output:
        zarr=directory(mcfg.out_dir / f'{params.wildcard_pattern}.zarr'),
    params:
        layer=lambda wildcards: mcfg.get_from_parameters(wildcards, 'counts', default='X'),
    localrule: True
    conda:
        get_env(config, 'scanpy')
    resources:
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb')
    script:
        '../scripts/collect.py'
