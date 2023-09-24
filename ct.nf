#!/usr/bin/env caastools

// I don't know if above env has to be set to nextflow

/* 
 * This code enables the new dsl of Nextflow. 
 */

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
single_alignment            : ${params.alignment}
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
discoveryOutput = "output_discovery"
resampleOutput = "output_resample"
bootstrapOutput = "ouptut_bootstrap"

single_alignment = file(params.alignment)
config = file(params.traitfile)

include { ct_discovery } from "${baseDir}/nf_modules/ct_discovery" addParams(OUTPUT: discoveryOutput, LABEL:"twocpus") 
//include { ct_resample } from "${baseDir}/ct" addParams(OUTPUT: resampleOutput, LABEL:"twocpus") 
//include { ct_bootstrap } from "${baseDir}/ct" addParams(OUTPUT: bootstrapOutput, LABEL:"twocpus") 

workflow {
	discovery_out = ct_discovery(single_alignment, config)
}



workflow.onComplete {
	println ( workflow.success ? "\nYay!\n" : "Oops .. something went wrong" )
}
