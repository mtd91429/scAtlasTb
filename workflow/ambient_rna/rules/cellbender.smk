def get_gpu_resource(resource_key, attempt):
    # stay on the GPU profile for retries: CellBender on CPU takes hours per batch
    return mcfg.get_resource(profile='gpu', resource_key=resource_key, attempt=attempt, attempt_to_cpu=attempt)


rule cellbender:
    """
    Run CellBender remove-background on the raw droplet matrix of one batch
    """
    input:
        batch=get_batch_file,
    output:
        h5=mcfg.out_dir / 'scatter' / params.wildcard_pattern / 'cellbender' / '{batch}' / 'cellbender.h5',
    params:
        args=lambda wildcards: get_methods(wildcards).get('cellbender', {}),
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
