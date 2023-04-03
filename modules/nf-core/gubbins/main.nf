process GUBBINS {
    tag "MAIN_GUBBINS"
    label 'process_high'

    conda (params.enable_conda ? 'gubbins' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.2--h79da9fb_0' :
        'quay.io/biocontainers/fastp:0.23.2--h79da9fb_0' }"

    input:
    path(contig_files)

    output:
    path('gubbins.final*.tre')             ,  emit: tree
    path('gubbins.node*.tre')              ,  emit: nodetree
    path('*predictions.gff')               ,  emit: predictions
    path('*statistics.csv')                ,  emit: stats
    path('*distribution.vcf')              ,  emit: vcf
    path('*.log')                          ,  emit: log
    path "versions.yml"                    ,  emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    """
    for file in *fa; do 
        echo -e "\$(basename \$file | cut -d. -f1)\t\$file"
    done > fasta_files.fofn

    generate_ska_alignment.py --threads ${task.cpus} --reference ${params.reference_fasta} --fasta fasta_files.fofn --out ${params.run_name}_align.fa

    run_gubbins.py --threads ${task.cpus} --prefix gubbins ${params.run_name}_align.fa

    mkdir -p gbb_output
    cp gubbins* gbb_output
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gubbins: \$(echo \$(run_gubbins.py --version 2>&1))
    END_VERSIONS
    """
}
