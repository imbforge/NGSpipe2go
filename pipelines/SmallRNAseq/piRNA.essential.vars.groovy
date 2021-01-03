//Pipeline generated with command line: ./imb-pip.pl --tasks-pip=1 --force
//By: ssayolsp At: Fr 17 Okt 2014 17:12:41 CEST
//
// REMEMBER TO CHANGE THESE ESSENTIAL VARS!!

/*COMMENT OUT THE SPECIES TO BE ANALYSED*/
// ZEBRAFISH
ESSENTIAL_PROJECT="/fsimb/groups/imb-kettinggr/adomingues/projects/test_pipeline"
ESSENTIAL_BOWTIE_REF="/fsimb/groups/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv10/index/bowtie/1.2.1.1/GRCz10"
ESSENTIAL_GENOME_REF="/fsimb/groups/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv10/Sequence/WholeGenomeFastaTopLevel/genome.fa"
ESSENTIAL_FEATURES="/fsimb/groups/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv10/Annotation/Genes/pipeline/transposons.bed"
ESSENTIAL_GENESGTF="/fsimb/groups/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv10/Annotation/Genes/pipeline/Danio_rerio_and_repeat_masker.GRCz10.80.chr.gtf"
ESSENTIAL_RRNA_BOWTIE_REF="/fsimb/groups/imb-kettinggr/adomingues/projects/test_pipeline/ref/drerio/rrna" // necessary for fastqscreen
ESSENTIAL_REPEAT_MASKER="/home/adomingu/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv10/Sequence/repEnrichNoSimpleLow"
ESSENTIAL_REPEAT_REF="/fsimb/groups/imb-kettinggr/adomingues/projects/test_pipeline/ref/drerio/RepEnrich"
REPENRICH_BED="FALSE" 
ESSENTIAL_THREADS=4
ESSENTIAL_SPECIES="zebrafish"   // necessary for miRDeep2, used to refer to UCSC

// adapater trimming
ESSENTIAL_READLENGTH=51      // actual read length in original raw data (incl. insert, UMIs, adapter)
ESSENTIAL_MINREADLENGTH=28 // minimal length of reads to be kept
MAX_LENGTH=38 // maximal length of reads to be kept for further analysis
ESSENTIAL_ADAPTER_SEQUENCE="-a AGATCGGAAGAGCACACGTCT -a AAGCAGTGGTATCAACGCAGAGT -a CTGTCTCTTATACACATCT -a TGGAATTCTCGGGTGCCAAGG -a CTCGTATGCCGTCTTCTGCTTG -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  // adapter to be trimmed off
ESSENTIAL_MINADAPTEROVERLAP=5 // minimal overlap of the read and the adapter

//UMI trimming
LEFT_UMI_LENGTH=4      // length of 5' UMI
RIGHT_UMI_LENGTH=4      // length of 3' UMI

// mapping and read counting
ESSENTIAL_SPECIES="zebrafish"   // necessary for miRDeep2, used to refer to UCSC
ESSENTIAL_SAMPLE_PREFIX="Sample_"
ESSENTIAL_FEATURETYPE="gene_type" //gencode uses gene_type; ensemble uses gene_biotype
ESSENTIAL_DB="danRer10"
ESSENTIAL_STRANDED="yes"
ESSENTIAL_PAIRED="no"
ESSENTIAL_MISMATCHES=2
ESSENTIAL_GENOME_SIZE=1400000000 // Taken from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1077594. Other sizes (see deppTools docs):  mm9: 2,150,570,000; hg19:2,451,960,000; dm3:121,400,000 and ce10:93,260,000


//global vars that will be reused in some global vars
PROJECT=ESSENTIAL_PROJECT
PROCESSED=PROJECT + "/rawdata_processed"
RESULTS=PROJECT + "/results"
LOGS=PROJECT + "/logs"
QC=RESULTS + "/qc"
REPORTS=RESULTS + "/reports"
MAPPED=RESULTS + "/mapped"
UNIQUE_MAPPED=MAPPED + "/unique"
UNMAPPED=MAPPED + "/unmapped"
FQ=MAPPED + "/bam2fq"
MULTI_MAPPED=MAPPED + "/multimapped"
TMP=PROJECT + "/tmp"
TRACKS=RESULTS + "/tracks"
