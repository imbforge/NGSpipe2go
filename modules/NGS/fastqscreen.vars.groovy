
//vars for task FastQC from catalog NGS, version 1
FASTQC_OUTDIR=QC + "/fastqc"
FASTQC_THREADS=4
//vars for task Fastqscreen
FASTQSCREEN_OUTDIR=QC + "/fastqscreen"
FASTQSCREEN_THREADS=4
//fastqscreen additional param e.g. subset or bowtie /bowtie 2 parameters
FASTQSCREEN_PARAM="--nohits --subset 100000"
//the fastqscreen_conf defines your references, with these we will create a fastqscreen conf script and then run the fastqscreen
//this could be e.g.
FASTQSCREEN_CONF="Human::BowtieIndex/genome,rRNA::mouse_all_rRNA.fasta"
