FASTQSCREEN_OUTDIR=QC + "/fastqscreen"
FASTQSCREEN_THREADS=Integer.toString(ESSENTIAL_THREADS)
//fastqscreen additional param e.g. subset or bowtie /bowtie 2 parameters
FASTQSCREEN_PARAM="--nohits --subset 100000"
//the fastqscreen_conf defines your references, with these we will create a fastqscreen conf script and then run the fastqscreen
//this could be e.g.
FASTQSCREEN_CONF="Human::BowtieIndex/genome,rRNA::mouse_all_rRNA.fasta"
