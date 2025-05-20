# uranus_nf_core_mini
A mini uranus mostly nf-core pipeline that incorporates mosdepth, cggpindel (not nf-core) and sentieon bwamem.

# Example Code to build the dnanexus applet
```
dx build --nextflow uranus_nf_core_mini --destination project-xxxx --extra-args '{"access":{"network": ["*"], "allProjects":"VIEW", "project":"UPLOAD"}}' --overwrite 
```
# Example Code to run this pipeline on DNAnexus:
```
dx run applet-xxxx -inextflow_run_opts='-profile docker' \
-inextflow_pipeline_params="--reads1='dx://project-xxxx:/Folder/file.fastq.gz' --reads2='dx://project-xxxx:/Folder/file.fastq.gz' 
--bwa_index='dx://project-xxxx:/Folder/file.bwa-index.tar.gz' 
--mosdepth_bed='dx://project-xxxx:/Folder/file.bed' 
--genomefa='dx://project-xxxx:/Folder/file.fasta.gz' 
--genomefa_fai='dx://project-xxxx:/Folder/file.fasta.gz.fai' 
--sample_name='file_name' 
--normal='dx://project-xxxx:/Folder/file.bam' 
--normal_index='dx://project-xxxx:/Folder/file.bam.bai' 
--genes='dx://project-xxxx:/Folder/file.bed' 
--unmatched='dx://project-xxxx:/Folder/file.gff3.gz' 
--simrep='dx://project-xxxx:/Folder/file.bed.gz' 
--simrep_index='dx://project-xxxx:/Folder/file.bed.gz.tbi' 
--filter='dx://project-xxxx:/Folder/file.lst'" 
```
