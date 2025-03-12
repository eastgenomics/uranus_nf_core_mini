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

    // Add error strategy for potential memory issues
    errorStrategy { task.exitStatus in [137,140] ? 'retry' : 'finish' }
    maxRetries 1

    input:
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
    publishDir '*', mode: 'copy'
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
    main:
    /*
    ------------------------------------------
    Step 1: Prepare input Channels
        ch_reads
        ch_bwa_index_archive
        ch_mosdepth_bed
        ch_fasta
    ------------------------------------------
    */
    ch_reads = Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .map { tuple -> [ [id: tuple[0]], tuple[1] ] }
    ch_bwa_index_archive = Channel.fromPath(params.bwa_index)    
    ch_mosdepth_bed = Channel.fromPath(params.mosdepth_bed)
    ch_fasta = Channel.value([[:], file(params.genomefa)])

    /*
    ------------------------------------------
    Step 2: Extract (untar) bwa index
    ------------------------------------------
    */

    EXTRACT_BWA_INDEX(ch_bwa_index_archive)
    ch_bwa_index = EXTRACT_BWA_INDEX.out.collect().map { files -> [[:], files] }

    /*
    ------------------------------------------
    Step 3: Alignment with BWA_MEM
    ------------------------------------------
    */
    BWA_MEM(
    ch_reads,
    ch_bwa_index, 
    ch_fasta, 
    false
    )   
    ch_bam_to_sort = BWA_MEM.out.bam.map { meta, bam ->
            def new_meta = meta + [id: "${meta.id}.sorted"]
            [ new_meta, bam ]
        }

    /*
    ------------------------------------------
    Step 4: Sort BAM Files with SAMTOOLS_SORT
    ------------------------------------------
    */
    SAMTOOLS_SORT(ch_bam_to_sort, ch_fasta)
    SAMTOOLS_SORT.out.bam.view { meta, bam -> 
        log.info "SAMTOOLS_SORT output: meta=${meta}, bam=${bam}"
        return [meta, bam]
    }

    /*
    ------------------------------------------
    Step 5: Index BAM Files with SAMTOOLS_INDEX
    ------------------------------------------
    */

    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)
    
    /*
    ------------------------------------------
    Step 6: Run Mosdepth
    ------------------------------------------
    */

    // Prepare inputs for mosdepth
    ch_bam_bai = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.bai)
    ch_mosdepth_input = ch_bam_bai.combine(ch_mosdepth_bed)
        .map { meta, bam, bai, bed -> 
            [ meta, bam, bai, bed ]
        }

    // Run MOSDEPTH
    MOSDEPTH(ch_mosdepth_input, ch_fasta)


    /*
    ------------------------------------------
    Step 7: Prepare inputs for cgppindel
    ------------------------------------------
    */
    ch_genomefa = Channel.fromPath(params.genomefa)
    ch_genomefa_fai = Channel.fromPath(params.genomefa_fai)
    ch_genes = Channel.fromPath(params.genes)
    ch_unmatched = Channel.fromPath(params.unmatched)
    ch_simrep = Channel.fromPath(params.simrep)
    ch_simrep_index = Channel.fromPath(params.simrep_index)
    ch_filter = Channel.fromPath(params.filter)



    // Using the BAM and BAI files generated from the FASTQ files as the tumor sample
    // But not in a single channel as was done for mosdepth
    ch_tumor_bam = SAMTOOLS_SORT.out.bam.map { meta, bam -> bam }
    ch_tumor_bai = SAMTOOLS_INDEX.out.bai.map { meta, bai -> bai }

    // Channels for the normal sample BAM and BAI files
    ch_normal_bam = Channel.fromPath(params.normal)
    ch_normal_bai = Channel.fromPath(params.normal_index)


    /*
    ------------------------------------------
    Step 8: run cgppindel
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


}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
