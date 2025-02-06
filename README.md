# uranus_nf_core_mini
A mini uranus mostly nf-core pipeline that incorporates mosdepth, cggpindel (not nf-core) and sentieon bwamem.


# Example Code to run this pipeline:
nextflow run main.nf -resume \
--reads 'XXX.fastq.gz' \
--bwa_index 'XXX.bwa-index.tar' \
--mosdepth_bed 'XXX.bed' \
--outdir '/output' \
--normal 'XXX.bam' \
--normal_index 'XXX.bam.bai' \
--tumour 'XXX.bam' \
--tumour_index 'XXX.bam.bai' \
--genomefa 'XXX.fa' \
--genomefa_fai 'XXX.fa.fai' \
--genes 'XXX.bed' \
--unmatched 'XXX.gff3.gz' \
--simrep 'XXX.bed.gz' \
--simrep_index 'XXX.bed.gz.tbi' \
--filter 'XXX.lst' -profile docker,laptop

