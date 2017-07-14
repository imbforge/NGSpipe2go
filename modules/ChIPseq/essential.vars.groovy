ESSENTIAL_PROJECT="/project"
ESSENTIAL_BOWTIE_REF="/annotation/mm9/Sequence/BowtieIndex/genome"
ESSENTIAL_BOWTIE_GENOME="/annotation/mm9/Sequence/BowtieIndex/genome.fa"
ESSENTIAL_CHR_SIZES=ESSENTIAL_PROJECT + "/mm9.chrom.sizes"
ESSENTIAL_SAMPLE_PREFIX=""
ESSENTIAL_BSGENOME="BSgenome.Mmusculus.UCSC.mm9"
ESSENTIAL_ANNODB="org.Mm.eg.db"   // it's not mandatory, but dont forget it if you wanna annotate your peaks!
ESSENTIAL_FRAGLEN=200
ESSENTIAL_READLEN=45
ESSENTIAL_MACS2_GSIZE="mm"
ESSENTIAL_THREADS=4
ESSENTIAL_FRAGMENT_USAGE="no" // essential variable which will tell bamCoverage to reconstitute fragments to create bigwig tracks
			      //set to yes in case of paired end DNA/ChIP sequencing usage
ESSENTIAL_PAIRED="no"	     //relevant for MACS2 if set it will used the PE mode in MACS2 peakcalling

ESSENTIAL_DUP="auto" // relevant for MACS2 it instructs macs2 to use its auto function to keep duplicate marked reads
//global vars that will be reused in some global vars
PROJECT=ESSENTIAL_PROJECT
LOGS=PROJECT + "/logs"
MAPPED=PROJECT + "/mapped"
QC=PROJECT + "/qc"
REPORTS=PROJECT + "/reports"
RESULTS=PROJECT + "/results"
TMP=PROJECT + "/tmp"
TRACKS=PROJECT + "/tracks"
