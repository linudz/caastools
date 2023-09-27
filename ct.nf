#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/* /*
 * CAASTools test pipe
 * @authors
 * Fabio Bartieri <fabio.barteri@upf.edu>
 * @collaborators
 * Miguel Ramon <miguel.ramon@upf.edu
 *
 */

/*
 *
 */

version = "0.0.1"
// this prevents a warning of undefined parameter
params.help             = false

// this prints the input parameters
log.info """
BIOCORE@CRG - N F TESTPIPE  ~  version ${version}
=============================================
alignment                   : ${params.alignment}
fmt                         : ${params.ali_format}
config                      : ${params.traitfile}
"""

// this prints the help in case you use --help parameter in the command line and it stops the pipeline
// this has and should be personalized or otherwise removed
if (params.help) {
    log.info 'Hi!'
    log.info 'Enjoy!'
    log.info '\n'
    exit 1
}

/*
 * Defining the output folders
 */

align_tuple = Channel
                .fromPath(params.alignment)
                .map { file -> tuple(file.baseName, file) }
config = file(params.traitfile)

include { ct_discovery } from "${baseDir}/subworkflows/ct_discovery" addParams(ALIGN_TUPLE: align_tuple, LABEL:"twocpus")

workflow {
	discovery_out = ct_discovery(align_tuple, config)
}



workflow.onComplete {
	println ( workflow.success ? "\nYay!\n" : "Oops .. something went wrong" )
}