//Pipeline generated with command line: ./imb-pip.pl --tasks-pip=1 --force
//By: ssayolsp At: Fr 17 Okt 2014 17:12:41 CEST
//
// REMEMBER TO CHANGE THESE ESSENTIAL VARS!!

// // ZEBRAFISH
// ESSENTIAL_PROJECT="/local/scratch1/imb-kettinggr/adomingues/projects/bpipe_small_rna"
// ESSENTIAL_BOWTIE_REF="/home/adomingu/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv9/Sequence/BowtieIndexChr/chr"
// ESSENTIAL_GENOME_REF="/home/adomingu/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv9/Sequence/chr_sequences/chr.clean.fa"
// // vars for piRNA analyis
// ESSENTIAL_FEATURES="/fsimb/groups/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv9/Annotation/pipeline/features.bed"
// // vars for mirDeep2
// ESSENTIAL_MATURE_MIRNA="~/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv9/Annotation/SmallRNA/mature.noSpaces.fa"
// ESSENTIAL_HAIRPIN_MIRNA="~/imb-git addkettinggr/genomes/Danio_rerio/Ensembl/Zv9/Annotation/SmallRNA/hairpin.noSpaces.fa"
// // read size to keep (added 8 bp to account for barcodes)
// // Zebrafish
// MIN_LENGTH=21
// MAX_LENGTH=51
// ESSENTIAL_STRANDED="yes"
// ESSENTIAL_PAIRED="no"

// C. ELEGANS
ESSENTIAL_PROJECT="/home/adomingu/imb-kettinggr/adomingues/projects/imb_ketting_2014_14_almeida_smallRNA_celegans"
ESSENTIAL_BOWTIE_REF="/home/adomingu/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/BowtieIndex/genome"
ESSENTIAL_GENOME_REF="/home/adomingu/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/BowtieIndex/genome.fa"
// vars for piRNA analyis
ESSENTIAL_GENESGTF="/home/adomingu/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.gtf"
ESSENTIAL_FEATURES="/home/adomingu/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.gtf"
ESSENTIAL_BIOTYPES_TABLE="/home/adomingu/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/wormbaseGeneID2biotype.txt"
// read size to keep (added 8 bp to account for barcodes)
MIN_LENGTH=26
MAX_LENGTH=38
ESSENTIAL_STRANDED="yes"
ESSENTIAL_PAIRED="no"
ESSENTIAL_THREADS=4

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

