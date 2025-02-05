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






process runCgpPindel {
    tag "$tumour"

    container 'quay.io/wtsicgp/cgppindel:3.9.0'
    containerOptions "--volume ${params.project_dir}:/data_two"

    // Set resource limits
    cpus 7
    memory 14.GB
    time '2h'  // Keep the time limit at 24 hours, adjust if needed

    // Add error strategy for potential memory issues
    //errorStrategy { task.exitStatus in [137,140] ? 'retry' : 'finish' }
    //maxRetries 1

    input:
    //val cgppindel_id
    path genomefa
    path genomefa_fai
    path simrep
    path simrep_index
    path genes
    path unmatched
    path filter
    path tumour
    path tumour_index
    path normal
    path normal_index
    val assembly
    val seqtype

    output:
    path "*.flagged.vcf.gz", emit: flagged_vcf
    path "*.flagged.vcf.gz.tbi", emit: vcf_index
    //path "output_vcf_with_vaf/*.af.vcf.gz", emit: vcf_with_vaf
    //path "output_log/logs.tar.gz", emit: logs
    //path "versions.yml", emit: versions
    publishDir '/home/raymondmiles/Desktop/Software_Development/software_module_master/eastgenomics-mini_uranus/', mode: 'copy'
    script:
    """
    mkdir -p output_vcf output_vcf_with_vaf output_log temp_logs
    ls -l ${genomefa}*
    echo "Running in directory: \$(pwd)"
    pindel.pl \\
    -reference ${genomefa} \\
    -simrep ${simrep} \\
    -genes ${genes} \\
    -unmatched ${unmatched} \\
    -assembly ${assembly} \\
    -species ${params.species} \\
    -seqtype ${seqtype} \\
    -filter ${filter} \\
    -tumour ${tumour} \\
    -normal ${normal} \\
    -outdir .

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

workflow MINI_URANUS {


    /*
    ------------------------------------------
    Step 1: Prepare Input Channels
    ------------------------------------------
    */


    /// take:
    /// ch_samplesheet // channel: samplesheet read in from --input
    main:
    //ch_versions = Channel.empty()
    ch_reads = Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .map { tuple -> [ [id: tuple[0]], tuple[1] ] }
        
    ch_bwa_index_archive = Channel.fromPath(params.bwa_index)
    mosdepth_bed = Channel.fromPath(params.mosdepth_bed) //mosdepth_bed.view{ "MOSDEPTH BED file: $it" }
    ch_fasta = Channel.value([[:], file(params.genomefa)])


    // Extract (untar) bwa index
    EXTRACT_BWA_INDEX(ch_bwa_index_archive)
    ch_bwa_index = EXTRACT_BWA_INDEX.out.collect().map { files -> [[:], files] }

    /*
    ------------------------------------------
    Step 2: Alignment with BWA MEM
    ------------------------------------------
    */
    BWA_MEM(
    ch_reads,
    ch_bwa_index, 
    ch_fasta, 
    false)

    ch_bam_to_sort = BWA_MEM.out.bam.map { meta, bam ->
            def new_meta = meta + [id: "${meta.id}.sorted"]
            [ new_meta, bam ]
        }
    /*
    ------------------------------------------
    Step 3: Sort BAM Files with SAMTOOLS_SORT
    ------------------------------------------
    */
    SAMTOOLS_SORT(ch_bam_to_sort, ch_fasta)
    SAMTOOLS_SORT.out.bam.view { meta, bam -> 
        log.info "SAMTOOLS_SORT output: meta=${meta}, bam=${bam}"
        return [meta, bam]
    }

    /*
    ------------------------------------------
    Step 4: Index BAM Files with SAMTOOLS_INDEX
    ------------------------------------------
    */
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)
    SAMTOOLS_INDEX.out.bai.view { meta, bai -> 
        log.info "SAMTOOLS_INDEX output: meta=${meta}, bai=${bai}"
        return [meta, bai]
    }


    /*
    ------------------------------------------
    Step 5: Run MOSDEPTH
    ------------------------------------------
    */
    // Prepare input for MOSDEPTH

    ch_bam_bai = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.bai)
    ch_bam_bai.view { meta, bam, bai ->
        log.info "Joined BAM and BAI: meta=${meta}, bam=${bam}, bai=${bai}"
        return [meta, bam, bai]
    }
    ch_mosdepth_input = ch_bam_bai.combine(mosdepth_bed)
    .map { meta, bam, bai, bed -> 
       log.info "MOSDEPTH input: meta=${meta}, bam=${bam}, bai=${bai}, bed=${bed}"
        [ meta, bam, bai, bed ]
    }

    // Run MOSDEPTH
    //MOSDEPTH(ch_mosdepth_input, ch_fasta)

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

        // Prepare inputs for runCgpPindel
        
/*
    ------------------------------------------
    Step 6: Prepare Inputs for runCgpPindel
    ------------------------------------------
    */

    // Channels for additional inputs
    ch_genomefa = Channel.fromPath(params.genomefa)
    ch_genomefa_fai = Channel.fromPath(params.genomefa_fai)
    ch_genes = Channel.fromPath(params.genes)
    ch_unmatched = Channel.fromPath(params.unmatched)
    ch_simrep = Channel.fromPath(params.simrep)
    ch_simrep_index = Channel.fromPath(params.simrep_index)
    ch_filter = Channel.fromPath(params.filter)

    // Using the BAM and BAI files generated from the FASTQ files as the tumor sample
ch_tumor_bam = SAMTOOLS_SORT.out.bam.map { meta, bam -> bam }
ch_tumor_bai = SAMTOOLS_INDEX.out.bai.map { meta, bai -> bai }

    // Channels for the normal sample BAM and BAI files
    ch_normal_bam = Channel.fromPath(params.normal)
    ch_normal_bai = Channel.fromPath(params.normal_index)

    /*
    ------------------------------------------
    Step 7: Run cgpPindel
    ------------------------------------------
    */

    runCgpPindel(
        ch_genomefa,
        ch_genomefa_fai,
        ch_simrep,
        ch_simrep_index,
        ch_genes,
        ch_unmatched,
        ch_filter,
        ch_tumor_bam,
        ch_tumor_bai,
        ch_normal_bam,
        ch_normal_bai,
        params.assembly,
        params.seqtype
    )


    emit:
        //bam = BWA_MEM.out.bam
        ///bam = SAMTOOLS_SORT.out.bam
        //bai = SAMTOOLS_INDEX.out.bai
        //mosdepth_global = MOSDEPTH.out.global_txt
        //mosdepth_regions = MOSDEPTH.out.regions_txt
        //mosdepth_thresholds = MOSDEPTH.out.thresholds_bed

        /// from cgppindel docker
        cgppindel_vcf = runCgpPindel.out.flagged_vcf
        cgppindel_vcf_index = runCgpPindel.out.vcf_index
     //   cgppindel_vcf_with_vaf = runCgpPindel.out.vcf_with_vaf
      //  cgppindel_logs = runCgpPindel.out.logs
        ///versions = BWA_MEM.out.versions.mix(MOSDEPTH.out.versions).mix(runCgpPindel.out.versions)

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
