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
    echo "DEBUG: Checking input file..."
    tar -xvf $bwa_index_archive || { echo "ERROR: tar extraction failed!" >&2; exit 1; }
    echo "DEBUG: Listing extracted files..."
    pwd
    ls -lh
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
        ch_fasta
    ------------------------------------------
    */


     if (!params.reads || !params.bwa_index || !params.genomefa) {
    error "ERROR: Missing required parameters. Please provide --reads, --bwa_index, and --genomefa."
    }
    
    /*************************
    
    View Parameters

    *************************/
    log.info "params.reads:  ${params.reads}"
    log.info "params.bwa_index:  ${params.bwa_index}"
    log.info "params.genomefa:  ${params.genomefa}"



    //Define channels and view channels
    ch_reads = Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .map { tuple -> [ [id: tuple[0]], tuple[1] ] }
        .view { "ch_reads: ${it}" } // View channel content
    ch_bwa_index_archive = Channel
        .fromPath(params.bwa_index)
        .view { "ch_bwa_index_archive: ${it}" } // View channel content

    ch_fasta = Channel
        .fromPath(params.genomefa)//Channel.value([[:], file(params.genomefa)])
        .view { "ch_fasta: ${it}" } // View channel content



    /*
    ------------------------------------------
    Step 2: Extract (untar) bwa index
    ------------------------------------------
    */
    log.info "starting extraction"

    EXTRACT_BWA_INDEX(ch_bwa_index_archive)
    ch_bwa_index = EXTRACT_BWA_INDEX.out.index_files
        .map { file -> tuple([:], file) }
        .view { "ch_bwa_index: ${it}" } // View channel content

    //ch_bwa_index = EXTRACT_BWA_INDEX.out.collect().map { files -> [[:], files] }
    
    log.info "extraction finished"

    /*
    ------------------------------------------
    Step 3: Alignment with BWA_MEM
    ------------------------------------------
    */
    log.info "BWA_MEM starting"

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
        .view { "ch_bam_to_sort: ${it}" } // View channel content

    log.info "BWA_MEM finished running"

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
    log.info "SAMTOOLS_SORT.out.bam:  ${SAMTOOLS_SORT.out.bam}"

    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)
    SAMTOOLS_INDEX.out.bai.view { meta, bai -> 
        log.info "SAMTOOLS_INDEX output: meta=${meta}, bai=${bai}"
        return [meta, bai]
    }
    emit:
        bam = SAMTOOLS_SORT.out.bam
        bai = SAMTOOLS_INDEX.out.bai

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
