process ROARY {
    tag "MAIN_ROARY"
    label 'process_high'

    conda (params.enable_conda ? 'roary' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.2--h79da9fb_0' :
        'quay.io/biocontainers/fastp:0.23.2--h79da9fb_0' }"

    input:
    path(gff_files)

    output:
    path('*')                          ,  emit: vf
    path('summary_statistics.txt')     , emit: summary
    path("versions.yml")                            , emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    """
    roary -p ${task.cpus} -r -e -n -qc ${gff_files} 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        roary: \$(echo \$(roary --version 2>&1))
    END_VERSIONS
    """
}
