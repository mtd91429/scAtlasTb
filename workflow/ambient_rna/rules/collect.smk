def collect_inputs(wildcards):
    methods = get_methods(wildcards)
    files = {
        'anndata': mcfg.get_input_file(**wildcards),
        'batches': expand_per_batch(f'{get_batch_dir(wildcards)}/{{batch}}.tsv', wildcards),
    }
    if 'cellbender' in methods:
        files['cellbender'] = expand_per_batch(rules.cellbender.output.h5, wildcards)
    return files


rule collect:
    input:
        unpack(collect_inputs)
    output:
        zarr=directory(mcfg.out_dir / f'{params.wildcard_pattern}.zarr'),
        summary=mcfg.image_dir / params.wildcard_pattern / 'summary.tsv',
    params:
        counts=lambda wildcards: mcfg.get_from_parameters(wildcards, 'counts', default='X'),
        batch_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'batch_key'),
        barcode_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'barcode'),
        gene_id_column=lambda wildcards: mcfg.get_from_parameters(wildcards, 'gene_id_column'),
    conda:
        get_env(config, 'scanpy')
    resources:
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='cpu', resource_key='mem_mb', attempt=attempt),
    script:
        '../scripts/collect.py'
