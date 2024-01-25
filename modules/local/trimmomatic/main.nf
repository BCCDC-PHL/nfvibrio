process TRIMMOMATIC {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::trimmomatic=0.39' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.2--h79da9fb_0' :
        'quay.io/biocontainers/fastp:0.23.2--h79da9fb_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.paired.fastq.gz') ,  emit: reads
    tuple val(meta), path('*.log')            , emit: log
    tuple val(meta), path('*.txt')            , emit: summary
    path "versions.yml"                       , emit: versions
    tuple val(meta), path('*.unpaired.fastq.gz')  , optional:true, emit: reads_fail

    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    [ ! -f  ${prefix}_1.fastq.gz ] && ln -sf ${reads[0]} ${prefix}_1.fastq.gz
    [ ! -f  ${prefix}_2.fastq.gz ] && ln -sf ${reads[1]} ${prefix}_2.fastq.gz

    trimmomatic PE \
    -threads ${task.cpus} \
    -trimlog ${prefix}_trimmomatic.log \
    -summary ${prefix}_summary.txt \
    ${prefix}_1.fastq.gz \
    ${prefix}_2.fastq.gz \
    ${prefix}_1.paired.fastq.gz \
    ${prefix}_1.unpaired.fastq.gz \
    ${prefix}_2.paired.fastq.gz \
    ${prefix}_2.unpaired.fastq.gz \
    SLIDINGWINDOW:4:25 MINLEN:25 TRAILING:25

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimmomatic: \$(echo \$(trimmomatic -version 2>&1))
    END_VERSIONS
    """
}
