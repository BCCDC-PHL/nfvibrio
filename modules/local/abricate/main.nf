process ABRICATE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? 'abricate' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.2--h79da9fb_0' :
        'quay.io/biocontainers/fastp:0.23.2--h79da9fb_0' }"

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path('*.vf.tsv')         ,  emit: vf
    // tuple val(meta), path('*.log')            , emit: log
    path "versions.yml"                       , emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    abricate --db vfdb ${contigs} > ${prefix}.vf.tsv 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abricate: \$(echo \$(abricate --version 2>&1) | sed 's/abricate //')
    END_VERSIONS
    """
}
