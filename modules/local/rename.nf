process RENAME {
    tag "$meta.id"

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
