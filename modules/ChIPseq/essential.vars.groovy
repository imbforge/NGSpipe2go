ESSENTIAL_PROJECT="/project"
ESSENTIAL_BOWTIE_REF="/annotation/mm9/Sequence/BowtieIndex/genome"
ESSENTIAL_BOWTIE_GENOME="/annotation/mm9/Sequence/BowtieIndex/genome.fa"
ESSENTIAL_CHROMSIZES="/annotation/mm9/mm9.chrom.sizes"  // chromosome sizes file
ESSENTIAL_SAMPLE_PREFIX=""
ESSENTIAL_BLACKLIST="/annotation/mm9/consensusBlacklist.bed" 
// Note that blacklist regions are applicable to human (hg38 and hg19), mouse (mm10 and mm9), worm (ce10) and fly (dm3) genomes.
ESSENTIAL_BSGENOME="BSgenome.Mmusculus.UCSC.mm9"
ESSENTIAL_TXDB="TxDb.Mmusculus.UCSC.mm9.knownGene" // it's not mandatory, but dont forget it if you wanna annotate your peaks
ESSENTIAL_ANNODB="org.Mm.eg.db"                    // it's not mandatory, but dont forget it if you wanna annotate your peaks
ESSENTIAL_FRAGLEN=200
ESSENTIAL_READLEN=45
ESSENTIAL_MACS2_GSIZE="mm"
ESSENTIAL_THREADS=4
ESSENTIAL_DB="mm9"             // UCSC assembly version for GREAT analysis. Note that GREAT only supports human (UCSC hg19 and UCSC hg18), 
                               // mouse (UCSC mm9) and zebrafish (UCSC danRer7) genomes.    
ESSENTIAL_FRAGMENT_USAGE="no" // essential variable which will tell bamCoverage to reconstitute fragments to create bigwig tracks
			      //set to yes in case of paired end DNA/ChIP sequencing usage
ESSENTIAL_PAIRED="no"	      //relevant for MACS2 if set it will used the PE mode in MACS2 peakcalling
ESSENTIAL_STRANDED="no"       // strandness: no|yes|reverse

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
