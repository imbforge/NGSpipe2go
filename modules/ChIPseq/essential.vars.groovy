ESSENTIAL_PROJECT="/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge/chipseq_module_test" // put project folder here
ESSENTIAL_BOWTIE_REF="/fsimb/groups/imb-bioinfocf/common-data/igenomes_reference/Mus_musculus/UCSC/mm9/Sequence/BowtieIndex/genome"
ESSENTIAL_CHR_SIZES=ESSENTIAL_PROJECT + "/mm9.chrom.sizes"
ESSENTIAL_SAMPLE_PREFIX="Sample_imb_gcf_2014_07_"
ESSENTIAL_BSGENOME="BSgenome.Mmusculus.UCSC.mm9"
ESSENTIAL_ANNODB="org.Mm.eg.db"   // it's not mandatory, but dont forget it if you wanna annotate your peaks!
ESSENTIAL_FRAGLEN=200
ESSENTIAL_READLEN=50
ESSENTIAL_MACS2_GSIZE="mm"
ESSENTIAL_THREADS=4

//global vars that will be reused in some global vars
PROJECT=ESSENTIAL_PROJECT
LOGS=PROJECT + "/logs"
MAPPED=PROJECT + "/mapped"
QC=PROJECT + "/qc"
REPORTS=PROJECT + "/reports"
RESULTS=PROJECT + "/results"
TMP=PROJECT + "/tmp"
TRACKS=PROJECT + "/tracks"
