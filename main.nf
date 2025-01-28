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

    //workflow pindel_pl ()
//

        // Prepare inputs for pindel_pl


    log.info "input prep"

            //pindel_pl(
            //take:
//
            //main:
            ////unmatched: params.unmatched,
            //genomefa: params.genomefa, //genome_ch,
            //tumour: params.tumour,
            //normal: params.normal,
            //species: params.species,
            //assembly: params.assembly,
            //outdir: params.outdir,
            //filter: file(params.filter),
            //genes: tuple(file(params.genes), file("${params.genes}.tbi")),
            //unmatched: tuple(file(params.unmatched), file("${params.unmatched}.tbi")),
            //simrep: tuple(file(params.simrep), file("${params.simrep}.tbi")),
            //seqtype: params.seqtype
        //)//

    //    // Run pindel_pl workflow
// pindel_pl()

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
 ///   pindel_pl(
  ///      genomefa: params.genomefa,
  ///      tumour: params.tumour,
  ///      normal: params.normal,
  ///      species: params.species,
  ///      assembly: params.assembly,
  ///      outdir: params.outdir,
  ///      filter: file(params.filter),
  ///      genes: tuple(file(params.genes), file("${params.genes}.tbi")),
  ///      unmatched: tuple(file(params.unmatched), file("${params.unmatched}.tbi")),
  ///      simrep: tuple(file(params.simrep), file("${params.simrep}.tbi")),
  ///      seqtype: params.seqtype
  ///  )
    //
    
    //PINDEL_PL(channel.value(params.outdir), channel.value(params.genomefa), channel.value(params.tumour), channel.value(params.normal), channel.value(params.simrep), channel.value(params.filter), channel.value(params.genes), channel.value(params.unmatched))

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
