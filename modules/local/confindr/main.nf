process CONFINDR {
    tag "MAIN_CONFINDR"
    label 'process_ultra'

    conda (params.enable_conda ? 'confindr' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.2--h79da9fb_0' :
        'quay.io/biocontainers/fastp:0.23.2--h79da9fb_0' }"

    input:
    path(fastq_files)

    output:
    path('confindr_out')                         , emit: vf
    path('confindr_out/*_log.txt')               , emit: log
    path('confindr_out/std.txt')                 , emit: std
    path('confindr_out/*_report.csv')            , emit: report
    path "versions.yml"                          , emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ''
    """
    mkdir -p fastqs &&
    mv ${fastq_files} fastqs &&
    confindr.py -t ${task.cpus} -d ${params.confindr_db} -i fastqs -o confindr_out > std.txt 2>&1 &&
    mv std.txt confindr_out


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        confindr: \$(echo \$(confindr.py --version 2>&1| sed 's/ConFindr //')
    END_VERSIONS
    """
}
