FASTQSCREEN_OUTDIR=QC + "/fastqscreen"
FASTQSCREEN_THREADS=ESSENTIAL_THREADS
//fastqscreen additiona param e.g. subset or bowtie /bowtie 2 parameters
FASTQSCREEN_PARAM="--nohits --subset 100000"
//the fastqscreen_conf defines your references, with these we will create a fastqscreen conf script and then run the fastqscreen
FASTQSCREEN_CONF="Human::/fsimb/groups/imb-bioinfocf/common-data/igenomes_reference/Mus_musculus/Ensembl/GRCm38/Sequence/BowtieIndex/genome,rRNA::/fsimb/groups/imb-bioinfocf/common-data/rrna/mouse/mouse_all_rRNA.fasta"
