process PLOT_VF_HEATMAP {
    tag "VF_HEATMAP"
    label 'process_single'

    errorStrategy 'ignore'

    conda "${projectDir}/environments/plot.yaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.2--h79da9fb_0' :
        'quay.io/biocontainers/fastp:0.23.2--h79da9fb_0' }"

    input:
    path(vf_outputs)

    output:
    path("vf_heatmap.html")         ,  emit: main
    path "versions.yml"             ,  emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p virulence_factors
    cp -L ${vf_outputs} virulence_factors

    run_notebook.py \
    --notebook ${projectDir}/bin/make_vf_heatmap.ipynb \
    --output vf_heatmap.html \
    virulence_factors

    """
}
