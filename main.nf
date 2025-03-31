#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    eastgenomics/mini_uranus_just_bwamem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/eastgenomics/mini_uranus_just_bwamem
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MINI_URANUS  } from './workflows/mini_uranus'
include {pindel_pl} from './workflows/cgppindel_main'


/// include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_mini_uranus_just_bwamem_pipeline'
/// include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_mini_uranus_just_bwamem_pipeline'
/*

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow EASTGENOMICS_MINI_URANUS {

    take:
  ///  samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    MINI_URANUS() ///(/// samplesheet)

    log.info "starting input prep"
    // Extract BAM files from MINI_URANUS_JUST_BWAMEM output
    //bam_files = MINI_URANUS_JUST_BWAMEM.out.bam
    //bai_files = MINI_URANUS_JUST_BWAMEM.out.bai


    log.info "input prep"


}



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
///    PIPELINE_INITIALISATION (
///        params.version,
///        params.validate_params,
///        params.monochrome_logs,
///        args,
///        params.outdir,
///        params.input
///    )

    // Wrap parameters into channels



    //EASTGENOMICS_MINI_URANUS_JUST_BWAMEM()


    //
    // WORKFLOW: Run main workflow
    //
    EASTGENOMICS_MINI_URANUS (
        //PIPELINE_INITIALISATION.out.samplesheet
    )
    // SUBWORKFLOW: Run completion tasks
    //
///    PIPELINE_COMPLETION (
///        params.outdir,
///        params.monochrome_logs,
///        
///        
///    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
