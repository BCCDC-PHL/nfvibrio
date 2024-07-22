process RENAME {
    tag "$meta.id"

    publishDir "${params.outdir}/contigs/", pattern: "*${extension}", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(file)
	val extension

    output:
    tuple val(meta), path("${meta.id}${extension}")    

    script: 
    """
	cp -r ${file} ${meta.id}${extension}
    """
}
