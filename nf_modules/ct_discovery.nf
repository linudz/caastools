/*
*  discovery module
*/

params.CONTAINER = "miralnso/caastools-barebones:latest"
params.OUTPUT = "discovery_output"

process discovery {
    
    publishDir(params.OUTPUT, mode: 'copy')
    tag { "${alignment}" }
    //label (params.LABEL)

    input:
    path alignment
    path config

    output:
    path "${alignment}.output", emit: alignment_out

    script:
    """
    ct discovery -a ${alignment} -t ${config} -o ${alignment}.output --fmt ${params.ali_format}
    """
}

workflow ct_discovery {

    take: 
        alignment
        traitfile
    main:
        discovery(alignment, traitfile)
    emit:
        disc_out = discovery.out.alignment_out 

}