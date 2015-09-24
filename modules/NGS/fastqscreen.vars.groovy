
//vars for task FastQC from catalog NGS, version 1
FASTQC_OUTDIR=QC + "/fastqc"
FASTQC_THREADS=4
//vars for task Fastqscreen
FASTQSCREEN_OUTDIR=QC + "/fastqscreen"
FASTQSCREEN_THREADS=4
//fastqscreen additiona param e.g. subset or bowtie /bowtie 2 parameters
FASTQSCREEN_PARAM="--nohits --subset 100000"
//the fastqscreen_conf defines your references, with these we will create a fastqscreen conf script and then run the fastqscreen
FASTQSCREEN_CONF="Human::/fsimb/groups/imb-bioinfocf/common-data/igenomes_reference/Mus_musculus/Ensembl/GRCm38/Sequence/BowtieIndex/genome,rRNA::/fsimb/groups/imb-bioinfocf/common-data/rrna/mouse/mouse_all_rRNA.fasta"
