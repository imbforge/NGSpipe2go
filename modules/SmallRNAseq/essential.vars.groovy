//Pipeline generated with command line: ./imb-pip.pl --tasks-pip=1 --force
//By: ssayolsp At: Fr 17 Okt 2014 17:12:41 CEST
//
// REMEMBER TO CHANGE THESE ESSENTIAL VARS!!

/*COMMENT OUT THE SPECIES TO BE ANALYSED*/
// C. ELEGANS
ESSENTIAL_PROJECT="/fsimb/groups/imb-kettinggr/adomingues/projects/imb_ketting_2017_10_sensor_and_mutant"
ESSENTIAL_BOWTIE_REF="/home/adomingu/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/BowtieIndexSensor/genome_sensor"
ESSENTIAL_GENOME_REF="/home/adomingu/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/BowtieIndexSensor/genome_sensor.fa"
// vars for smRNA analysis
ESSENTIAL_GENES="/fsimb/groups/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.gtf"
ESSENTIAL_FEATURES="/fsimb/groups/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.gtf"
ESSENTIAL_TRANSPOSONS="/fsimb/groups/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/c_elegans.WS235.transposons.gtf"
ESSENTIAL_21U_REF="/fsimb/groups/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/biotypes/piRNA.gtf"
ESSENTIAL_22G_REF="/fsimb/groups/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/proteincoding_pseudogenes_lincRNA_transposons.pipeline.gtf"
ESSENTIAL_26G_REF="/fsimb/groups/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/proteincoding_pseudogenes_lincRNA.pipeline.gtf"
ESSENTIAL_REPEAT_MASKER="/fsimb/groups/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/repeatMasker/c_elegans.Wb235.ce11.repeatMasker.noSimpleLow.bed"
ESSENTIAL_REPEAT_REF="/fsimb/groups/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/repEnrichNoSimpleLow"
ESSENTIAL_BIOTYPES_TABLE="/fsimb/groups/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/wormbaseGeneID2biotype.txt"
NORMALIZATION_TO_NONSTRUCT="yes"
ESSENTIAL_ORG=""           // UCSC organism
ESSENTIAL_DB="" 
REPENRICH_BED="TRUE" // the annotation was prepared as bed file.
// read size to keep (added 8 bp to account for barcodes)
MIN_LENGTH=26
MAX_LENGTH=38
ESSENTIAL_STRANDED="reverse"
ESSENTIAL_PAIRED="no"
ESSENTIAL_THREADS=4
ESSENTIAL_MISMATCHES=0


//global vars that will be reused in some global vars
PROJECT=ESSENTIAL_PROJECT
LOGS=PROJECT + "/logs"
QC=PROJECT + "/data/qc"
REPORTS=PROJECT + "/reports"
RESULTS=PROJECT + "/results"
MAPPED=RESULTS + "/mapped"
UNIQUE_MAPPED=MAPPED + "/unique"
UNMAPPED=MAPPED + "/unmapped"
FQ=MAPPED + "/bam2fq"
MULTI_MAPPED=MAPPED + "/multimapped"
TMP=PROJECT + "/tmp"
TRACKS=MAPPED + "/tracks"
