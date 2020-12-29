//Pipeline generated with command line: ./imb-pip.pl --tasks-pip=1 --force
//By: ssayolsp At: Fr 17 Okt 2014 17:12:41 CEST
//
// REMEMBER TO CHANGE THESE ESSENTIAL VARS!!

/*COMMENT OUT THE SPECIES TO BE ANALYSED*/
// C. ELEGANS
ESSENTIAL_PROJECT="/fsimb/groups/imb-kettinggr/adomingues/projects/test_pipeline"
ESSENTIAL_BOWTIE_REF="/fsimb/groups/imb-kettinggr/adomingues/projects/test_pipeline/ref/cel_sensor"
ESSENTIAL_GENOME_REF="/fsimb/groups/imb-kettinggr/adomingues/projects/test_pipeline/ref/cel_sensor.fa"
// vars for smRNA analysis
ESSENTIAL_GENES="/fsimb/groups/imb-kettinggr/adomingues/projects/test_pipeline/ref/cel.gtf"
ESSENTIAL_FEATURES="/fsimb/groups/imb-kettinggr/adomingues/projects/test_pipeline/ref/cel.gtf"
ESSENTIAL_TRANSPOSONS="/fsimb/groups/imb-kettinggr/adomingues/projects/test_pipeline/ref/cel.transposons.gtf"
ESSENTIAL_21U_REF="/fsimb/groups/imb-kettinggr/adomingues/projects/test_pipeline/ref/piRNA.pipeline.gtf"
ESSENTIAL_22G_REF="/fsimb/groups/imb-kettinggr/adomingues/projects/test_pipeline/ref/proteincoding_pseudogenes_lincRNA_transposons.pipeline.gtf"
ESSENTIAL_26G_REF="/fsimb/groups/imb-kettinggr/adomingues/projects/test_pipeline/ref/proteincoding_pseudogenes_lincRNA.pipeline.gtf"
ESSENTIAL_MIRNA_REF="/fsimb/groups/imb-kettinggr/adomingues/projects/test_pipeline/ref/miRNA.pipeline.gtf"
ESSENTIAL_SENSOR_FILTER_REF="NGSpipe2go/testdata/siRNA/21U_sensor_full.bed"
ESSENTIAL_SENSOR_REF="NGSpipe2go/testdata/siRNA/21U_sensor_to_plot.bed"
ESSENTIAL_STUCTURAL_REF="/fsimb/groups/imb-kettinggr/adomingues/projects/test_pipeline/ref/structural.pipeline.gtf"
ESSENTIAL_ORG=""           // UCSC organism
ESSENTIAL_DB="" 

// adapater trimming
ESSENTIAL_READLENGTH=51      // actual read length in original raw data (incl. insert, UMIs, adapter)
ESSENTIAL_MINREADLENGTH=26 // minimal length of reads to be kept
MAX_LENGTH=38 // maximal length of reads to be kept for further analysis
ESSENTIAL_ADAPTER_SEQUENCE="-a AGATCGGAAGAGCACACGTCT -a AAGCAGTGGTATCAACGCAGAGT -a CTGTCTCTTATACACATCT -a TGGAATTCTCGGGTGCCAAGG -a CTCGTATGCCGTCTTCTGCTTG"  // adapter to be trimmed off
ESSENTIAL_MINADAPTEROVERLAP=5 // minimal overlap of the read and the adapter

//UMI trimming
LEFT_UMI_LENGTH=4      // length of 5' UMI
RIGHT_UMI_LENGTH=4      // length of 3' UMI

// mapping and read counting
ESSENTIAL_STRANDED="reverse"
ESSENTIAL_PAIRED="no"
ESSENTIAL_THREADS=4
ESSENTIAL_MISMATCHES=0


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
