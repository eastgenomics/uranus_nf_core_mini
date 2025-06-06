/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MERGE_BAM {
    tag "Merging BAM files"
    
    input:
    tuple val(meta1), path(bam1)
    tuple val(meta2), path(bam2)
    
    output:
    tuple val(meta1), path("merged.bam")
    
    script:
    """
    samtools merge merged.bam ${bam1} ${bam2}
    """
}


process EXTRACT_BWA_INDEX {
    errorStrategy 'retry'
    maxRetries 3

    input:
    path bwa_index_archive

    output:
    path '*', emit: index_files


    script:
    """
    tar -xvf $bwa_index_archive || { echo "ERROR: tar extraction failed!" >&2; exit 1; }
    """
}

process runCgpPindel {
    tag "$tumour"
    label 'process_high'

    container 'quay.io/wtsicgp/cgppindel:3.9.0'

    // Add error strategy for potential memory issues
    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    maxRetries 3

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
    script:
    """
    mkdir -p output_vcf output_vcf_with_vaf output_log temp_logs
    ls -l ${genomefa}*
    gunzip -c ${genomefa} > uncompressed.fa
    samtools faidx uncompressed.fa

    echo "Running in directory: \$(pwd)"
    pindel.pl \\
    -reference uncompressed.fa \\
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

include { BWA_MEM as BWA_MEM_LANE1 } from '../modules/nf-core/bwa/mem/main'
include { BWA_MEM as BWA_MEM_LANE2 } from '../modules/nf-core/bwa/mem/main'
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
        ch_fasta
    ------------------------------------------
    */
    
    log.info "params.reads1:  ${params.reads1}"
    log.info "params.reads2:  ${params.reads2}"
    log.info "params.bwa_index:  ${params.bwa_index}"
    log.info "params.genomefa:  ${params.genomefa}"


    ch_reads_lane1 = Channel
        .fromFilePairs(params.reads1, checkIfExists: true)
        .ifEmpty { error "No matching reads files found: ${params.reads1}" }
        .map { tuple -> [ [ id: tuple[0] ], tuple[1] ] }

    ch_reads_lane2 = Channel
        .fromFilePairs(params.reads2, checkIfExists: true)
        .ifEmpty { error "No matching reads files found: ${params.reads2}" }
        .map { tuple -> [ [ id: tuple[0] ], tuple[1] ] }

    ch_bwa_index_archive = Channel.fromPath(params.bwa_index)
        .ifEmpty { error "No matching bwa index files found: ${params.bwa_index}"}

    ch_fasta = Channel.fromPath(params.genomefa)
            .ifEmpty { error "Genome FASTA file not found: ${params.genomefa}" }
            .map { file -> [[:], file] }

    ch_mosdepth_bed = Channel.fromPath(params.mosdepth_bed)
            .ifEmpty { error "Mosdepth BED file not found: ${params.mosdepth_bed}" }

    /*
    ------------------------------------------
    Step 2: Extract (untar) bwa index
    ------------------------------------------
    */
    log.info "Starting BWA index extraction"
    EXTRACT_BWA_INDEX(ch_bwa_index_archive)
    ch_bwa_index = EXTRACT_BWA_INDEX.out.index_files.map { file -> [[:], file] }
    log.info "BWA index extraction finished"


    /*
    ------------------------------------------
    Step 3: Alignment with BWA_MEM
    ------------------------------------------
    */
    log.info "Starting BWA MEM process"

    // Run BWA_MEM separately on each lane
    BWA_MEM_LANE1(ch_reads_lane1, ch_bwa_index, ch_fasta, false)
    BWA_MEM_LANE2(ch_reads_lane2, ch_bwa_index, ch_fasta, false)

    /*
    ------------------------------------------
    Step 4: Join BAMs for merging
    ------------------------------------------
    */
    MERGE_BAM(BWA_MEM_LANE1.out.bam, BWA_MEM_LANE2.out.bam)
    ch_bam_to_sort = MERGE_BAM.out.map { meta, bam ->
            def new_meta = meta + [ id: "${meta.id}.sorted" ]
            [ new_meta, bam ]
        }

    log.info "Finished BWA MEM process"

    /*
    ------------------------------------------
    Step 5: Sort BAM Files with SAMTOOLS_SORT
    ------------------------------------------
    */
    log.info "Starting Samtools sorting and indexing"
    SAMTOOLS_SORT(ch_bam_to_sort, ch_fasta)
    SAMTOOLS_SORT.out.bam.view { meta, bam -> 
        log.info "SAMTOOLS_SORT output: meta=${meta}, bam=${bam}"
        return [meta, bam]
    }

    /*
    ------------------------------------------
    Step 6: Index BAM Files with SAMTOOLS_INDEX
    ------------------------------------------
    */
    log.info "SAMTOOLS_SORT.out.bam:  ${SAMTOOLS_SORT.out.bam}"

    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)
    SAMTOOLS_INDEX.out.bai.view { meta, bai -> 
        log.info "SAMTOOLS_INDEX output: meta=${meta}, bai=${bai}"
        return [meta, bai]
    }
    log.info "Finished Samtools sorting and indexing"

    /*
    ------------------------------------------
    Step 7: Run Mosdepth
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
    Step 8: Prepare inputs for cgppindel
    ------------------------------------------
    */
    ch_genomefa = Channel.fromPath(params.genomefa)
        .ifEmpty { error "Genome FASTA file not found: ${params.genomefa}" }
    ch_genomefa_fai = Channel.fromPath(params.genomefa_fai)
        .ifEmpty { error "Genome FASTA index file not found: ${params.genomefa_fai}" }
    ch_genes = Channel.fromPath(params.genes)
        .ifEmpty { error "Genes file not found: ${params.genes}" }
    ch_unmatched = Channel.fromPath(params.unmatched)
        .ifEmpty { error "Unmatched file not found: ${params.unmatched}" }
    ch_simrep = Channel.fromPath(params.simrep)
        .ifEmpty { error "Simrep file not found: ${params.simrep}" }
    ch_simrep_index = Channel.fromPath(params.simrep_index)
        .ifEmpty { error "Simrep index file not found: ${params.simrep_index}" }
    ch_filter = Channel.fromPath(params.filter)
        .ifEmpty { error "Filter file not found: ${params.filter}" }

    // Using the BAM and BAI files generated from the FASTQ files as the tumor sample
    // But not in a single channel as was done for mosdepth
    ch_tumor_bam = SAMTOOLS_SORT.out.bam.map { meta, bam -> bam }
    ch_tumor_bai = SAMTOOLS_INDEX.out.bai.map { meta, bai -> bai }

    // Channels for the normal sample BAM and BAI files
    ch_normal_bam = Channel.fromPath(params.normal)
        .ifEmpty { error "Normal BAM file not found: ${params.normal}" }
    ch_normal_bai = Channel.fromPath(params.normal_index)
        .ifEmpty { error "Normal BAM index file not found: ${params.normal_index}" }


    /*
    ------------------------------------------
    Step 9: run cgppindel
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
        bam = SAMTOOLS_SORT.out.bam
        bai = SAMTOOLS_INDEX.out.bai

}
workflow.onError {
    log.error "Pipeline execution stopped with error: ${workflow.errorMessage}"
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
