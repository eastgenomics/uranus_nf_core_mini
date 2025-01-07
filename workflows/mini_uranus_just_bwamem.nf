/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


process EXTRACT_BWA_INDEX {
    input:
    path bwa_index_archive

    output:
    path "*"

    script:
    """
    tar -xvf $bwa_index_archive
    """
}



///include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

include { BWA_MEM } from '../modules/nf-core/bwa/mem/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MINI_URANUS_JUST_BWAMEM {

    /// take:
    /// ch_samplesheet // channel: samplesheet read in from --input
    main:

    //ch_versions = Channel.empty()
    

        // Input channels
    ch_reads = Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .map { tuple -> [ [id: tuple[0]], tuple[1] ] }
    ch_bwa_index_archive = Channel.fromPath(params.bwa_index)

    //ch_index = Channel.value([[:], file(params.bwa_index)])
    ch_fasta = Channel.value([[:], file(params.fasta)])

    EXTRACT_BWA_INDEX(ch_bwa_index_archive)

    // Run BWA MEM
    BWA_MEM(ch_reads, EXTRACT_BWA_INDEX.out, ch_fasta, false)


    //
    // Collate and save software versions
    //
    /// softwareVersionsToYAML(ch_versions)
    ///     .collectFile(
    ///         storeDir: "${params.outdir}/pipeline_info",
    ///         name:  ''  + 'pipeline_software_' +  ''  + 'versions.yml',
    ///         sort: true,
    ///         newLine: true
    ///     ).set { ch_collated_versions }


    emit:
    //versions       = ch_versions                 // channel: [ path(versions.yml) ]
    bam = BWA_MEM.out.bam
    bam = BWA_MEM.out.bai


}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
