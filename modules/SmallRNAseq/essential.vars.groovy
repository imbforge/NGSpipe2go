//Pipeline generated with command line: ./imb-pip.pl --tasks-pip=1 --force
//By: ssayolsp At: Fr 17 Okt 2014 17:12:41 CEST
//
// REMEMBER TO CHANGE THESE ESSENTIAL VARS!!

/*COMMENT OUT THE SPECIES TO BE ANALYSED*/
// C. ELEGANS
/*ESSENTIAL_PROJECT="."
ESSENTIAL_BOWTIE_REF="/home/adomingu/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/BowtieIndexSensor/genome_sensor"
ESSENTIAL_GENOME_REF="/home/adomingu/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/BowtieIndexSensor/genome_sensor.fa"
// vars for smRNA analysis
ESSENTIAL_GENES="/fsimb/groups/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.gtf"
ESSENTIAL_FEATURES="/fsimb/groups/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.gtf"
ESSENTIAL_TRANSPOSONS="/fsimb/groups/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/c_elegans.WS235.transposons.gtf"
ESSENTIAL_21U_REF="/fsimb/groups/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/biotypes/piRNA.gtf"
ESSENTIAL_22G_REF="/fsimb/groups/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/proteincoding_pseudogenes_lincRNA_transposons.pipeline.gtf"
ESSENTIAL_26G_REF="/fsimb/groups/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/proteincoding_pseudogenes_lincRNA.pipeline.gtf"
ESSENTIAL_SENSOR_REF="/fsimb/groups/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Sensor/21U_sensor_to_plot.bed"
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
ESSENTIAL_MISMATCHES=0*/


// ZEBRAFISH
ESSENTIAL_PROJECT="/fsimb/groups/imb-kettinggr/adomingues/projects/imb_ketting_2018_01_redl_pgcs_total"
ESSENTIAL_STAR_REF="/fsimb/groups/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv10/index/star/2.5.2b"
ESSENTIAL_GENESGTF="/fsimb/groups/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv10/Annotation/Genes/Danio_rerio.GRCz10.80.chr.gtf"
ESSENTIAL_GENESBED="/fsimb/groups/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv10/Annotation/Genes/Danio_rerio.GRCz10.80.chr.bed"
ESSENTIAL_BOWTIE_REF="/fsimb/groups/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv10/index/bowtie/1.2.1.1/GRCz10"
ESSENTIAL_GENOME_REF="/fsimb/groups/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv10/Sequence/WholeGenomeFastaTopLevel/genome.fa"
// vars for piRNA analyis
ESSENTIAL_FEATURES="/fsimb/groups/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv10/Annotation/Genes/pipeline/transposons.bed"
ESSENTIAL_GENES="/fsimb/groups/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv10/Annotation/Genes/pipeline/Danio_rerio_and_repeat_masker.GRCz10.80.chr.gtf"
ESSENTIAL_REPEAT_MASKER="/home/adomingu/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv10/Sequence/repeatMasker/danRer10.noSimpleLow.fa.out"
ESSENTIAL_REPEAT_REF="/home/adomingu/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv10/Sequence/repEnrichNoSimpleLow"
ESSENTIAL_ADAPTER_SEQUENCE="AGATCGGAAGAGCACACGTCT"
REPENRICH_BED="FALSE" 
ESSENTIAL_THREADS=8
// read size to keep (added 8 bp to account for barcodes)
MIN_LENGTH=28
MAX_LENGTH=45
ESSENTIAL_STRANDED="yes"
ESSENTIAL_PAIRED="no"
ESSENTIAL_MISMATCHES=2
ESSENTIAL_DB="danRer10"
ESSENTIAL_GENOME_SIZE=1400000000 // Taken from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1077594. Other sizes (see deppTools docs):  mm9: 2,150,570,000; hg19:2,451,960,000; dm3:121,400,000 and ce10:93,260,000

/*// vars for mirDeep2
ESSENTIAL_MATURE_MIRNA="/fsimb/groups/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv9/Annotation/SmallRNA/mature.noSpaces.fa"
ESSENTIAL_HAIRPIN_MIRNA="~/imb-git addkettinggr/genomes/Danio_rerio/Ensembl/Zv9/Annotation/SmallRNA/hairpin.noSpaces.fa"*/

//global vars that will be reused in some global vars
PROJECT=ESSENTIAL_PROJECT
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
