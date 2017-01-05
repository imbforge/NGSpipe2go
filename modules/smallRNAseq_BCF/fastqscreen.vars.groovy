FASTQSCREEN_OUTDIR=QC + "/fastqscreen"
FASTQSCREEN_THREADS=Integer.toString(ESSENTIAL_THREADS)

//fastqscreen additional param e.g. subset or bowtie /bowtie 2 parameters
FASTQSCREEN_PARAM="--subset 100000 --aligner bowtie --bowtie \"-v2\""

//the fastqscreen_conf defines your references, with these we will create a fastqscreen conf script and then run the fastqscreen
//this could be e.g.
FASTQSCREEN_CONF=ESSENTIAL_SPECIES + "_genome::" + ESSENTIAL_BOWTIE_REF + "," + ESSENTIAL_SPECIES + "_rRNA::" + ESSENTIAL_RRNA_BOWTIE_REF
