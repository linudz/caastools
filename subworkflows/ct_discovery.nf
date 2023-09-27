/*
*  discovery module
*/

params.CONTAINER = "miralnso/caastools-barebones:latest"
params.OUTPUT = "results/${nextflow.timestamp}/discovery_output"

process discovery {
    // Define where to publish the output files.
    publishDir(params.OUTPUT, mode: 'copy')

    input:
    tuple val(alignmentID), file(alignmentFile)
    file traitfile

    output:
    tuple val(alignmentID), file("${alignmentID}.output")
    
    container params.CONTAINER  // Set the container here

    script:
    """
    ct discovery -a ${alignmentFile} -t ${params.traitfile} -o ${alignmentID}.output --fmt ${params.ali_format}
    """
}

workflow ct_discovery {
    take: 
        align_tuple
        traitfile
    main:
        discovery(align_tuple, traitfile)
    emit:
        disc_out = discovery.out
}