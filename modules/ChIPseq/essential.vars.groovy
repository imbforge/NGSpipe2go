ESSENTIAL_PROJECT="/fsimb/groups/imb-bioinfocf/projects/cfb_internal/test_projects/ChIP_seq_pe/"
ESSENTIAL_BOWTIE_REF="/fsimb/groups/imb-bioinfocf/common-data/igenomes_reference/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome"
ESSENTIAL_BOWTIE_GENOME="/fsimb/groups/imb-bioinfocf/common-data/igenomes_reference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
ESSENTIAL_CHR_SIZES=ESSENTIAL_PROJECT + "/hg19.chrom.sizes"
ESSENTIAL_SAMPLE_PREFIX="imb_roukos_2016_04_"
ESSENTIAL_BSGENOME="BSgenome.Hsapiens.UCSC.hg19"
ESSENTIAL_ANNODB="org.Hs.eg.db"   // it's not mandatory, but dont forget it if you wanna annotate your peaks!
ESSENTIAL_FRAGLEN=200
ESSENTIAL_READLEN=43
ESSENTIAL_MACS2_GSIZE="hs"
ESSENTIAL_THREADS=4
ESSENTIAL_FRAGMENT_USAGE="yes" // essential variable which will tell bamCoverage to reconstitute fragments to create bigwig tracks
			      //set to yes in case of paired end DNA/ChIP sequencing usage
ESSENTIAL_PAIRED="yes"	     //relevant for MACS2 if set it will used the PE mode in MACS2 peakcalling
//global vars that will be reused in some global vars
PROJECT=ESSENTIAL_PROJECT
LOGS=PROJECT + "/logs"
MAPPED=PROJECT + "/mapped"
QC=PROJECT + "/qc"
REPORTS=PROJECT + "/reports"
RESULTS=PROJECT + "/results"
TMP=PROJECT + "/tmp"
TRACKS=PROJECT + "/tracks"
