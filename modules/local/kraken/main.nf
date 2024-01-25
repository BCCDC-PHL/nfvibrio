process KRAKEN {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'kraken2 krakentools' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.2--h79da9fb_0' :
        'quay.io/biocontainers/fastp:0.23.2--h79da9fb_0' }"

    input:
    tuple val(meta), path(fastqs)

    output:
    tuple val(meta), path('*.report')            ,  emit: report
    tuple val(meta), path('*.out')               ,  emit: output
    tuple val(meta), path('*.kfilter.fastq.gz')  , emit: fastq
    path "versions.yml"                          , emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    kraken2 --confidence ${params.kraken_confidence} \
    --threads ${task.cpus} --db ${params.kraken_database} \
    --paired ${fastqs[0]} ${fastqs[1]} \
    --report ${prefix}.kraken.report > ${prefix}.kraken.out 
    
    extract_kraken_reads.py --fastq-output \
    -k ${prefix}.kraken.out -r ${prefix}.kraken.report \
    -1 ${fastqs[0]} -2 ${fastqs[1]} \
    -o ${prefix}_R1.kfilter.fastq -o2 ${prefix}_R2.kfilter.fastq \
    -t ${params.kraken_vibrio_id} 0 --include-children &&
    gzip -c ${prefix}_R1.kfilter.fastq > ${prefix}_R1.kfilter.fastq.gz &&
    gzip -c ${prefix}_R2.kfilter.fastq > ${prefix}_R2.kfilter.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(kraken2 --version | head -n1 | sed 's/Kraken version //g')
    END_VERSIONS
    """
}
