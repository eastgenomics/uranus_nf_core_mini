/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


process EXTRACT_BWA_INDEX {
    input:
    path bwa_index_archive

    output:
    path "*", emit: index_files


    script:
    """
    tar -xvf $bwa_index_archive
    """
}





///include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

include { BWA_MEM } from '../modules/nf-core/bwa/mem/main'
include { MOSDEPTH } from '../modules/nf-core/mosdepth/main'
include { SAMTOOLS_SORT } from '../modules/nf-core/samtools/sort/main'                                                                   
include { SAMTOOLS_INDEX } from '../modules/nf-core/samtools/index/main'                                                                 
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

    mosdepth_bed = Channel.fromPath(params.mosdepth_bed) //mosdepth_bed.view{ "MOSDEPTH BED file: $it" }

    ch_fasta = Channel.value([[:], file(params.fasta)])


    // Extract (untar) bwa index
    EXTRACT_BWA_INDEX(ch_bwa_index_archive)
    ch_bwa_index = EXTRACT_BWA_INDEX.out.collect().map { files -> [[:], files] }

    // Run BWA MEM
    BWA_MEM(ch_reads, ch_bwa_index, ch_fasta, false)
    ch_bam_to_sort = BWA_MEM.out.bam.map { meta, bam ->
            def new_meta = meta + [id: "${meta.id}.sorted"]
            [ new_meta, bam ]
        }

    // Sort BAM
    SAMTOOLS_SORT(ch_bam_to_sort, ch_fasta)
    //SAMTOOLS_SORT.out.bam.view { meta, bam -> 
    //    log.info "SAMTOOLS_SORT output: meta=${meta}, bam=${bam}"
    //    return [meta, bam]
    //}

    // Index BAM
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)
    //SAMTOOLS_INDEX.out.bai.view { meta, bai -> 
    //    log.info "SAMTOOLS_INDEX output: meta=${meta}, bai=${bai}"
    //    return [meta, bai]
    //}

    // Prepare input for MOSDEPTH
    ch_bam_bai = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.bai)
    //ch_bam_bai.view { meta, bam, bai ->
    //    log.info "Joined BAM and BAI: meta=${meta}, bam=${bam}, bai=${bai}"
    //    return [meta, bam, bai]
    //}

    // Prepare MOSDEPTH input
    ch_mosdepth_input = ch_bam_bai.combine(mosdepth_bed)
    .map { meta, bam, bai, bed -> 
        log.info "MOSDEPTH input: meta=${meta}, bam=${bam}, bai=${bai}, bed=${bed}"
        [ meta, bam, bai, bed ]
    }

    // Run MOSDEPTH
    MOSDEPTH(ch_mosdepth_input, ch_fasta)

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
    //bam = BWA_MEM.out.bam
    bam = SAMTOOLS_SORT.out.bam
    bai = SAMTOOLS_INDEX.out.bai
    mosdepth_global = MOSDEPTH.out.global_txt
    mosdepth_regions = MOSDEPTH.out.regions_txt
    mosdepth_thresholds = MOSDEPTH.out.thresholds_bed
    versions = BWA_MEM.out.versions.mix(MOSDEPTH.out.versions)

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
