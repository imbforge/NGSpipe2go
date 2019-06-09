//Pipeline generated with command line: ./imb-pip.pl --tasks-pip=1 --force
//By: ssayolsp At: Fr 17 Okt 2014 17:12:41 CEST
//
// REMEMBER TO CHANGE THESE ESSENTIAL VARS!!

/*COMMENT OUT THE SPECIES TO BE ANALYSED*/
// C. ELEGANS
/*ESSENTIAL_PROJECT="/project/"
ESSENTIAL_BOWTIE_REF="/home/adomingu/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/BowtieIndexSensor/genome_sensor"
ESSENTIAL_GENOME_REF="/home/adomingu/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/BowtieIndexSensor/genome_sensor.fa"
// vars for smRNA analysis
ESSENTIAL_GENESGTF="/fsimb/groups/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.gtf"
ESSENTIAL_RRNA_BOWTIE_REF="/data/rrna/celegans/celegans_all_rRNA" // necessary for fastqscreen

// needed to filter small RNA classes
ESSENTIAL_FEATURES="/fsimb/groups/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.gtf"
ESSENTIAL_TRANSPOSONS="/fsimb/groups/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/c_elegans.WS235.transposons.gtf"
ESSENTIAL_21U_REF="/fsimb/groups/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/biotypes/piRNA.gtf"
ESSENTIAL_22G_REF="/fsimb/groups/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/proteincoding_pseudogenes_lincRNA_transposons.pipeline.gtf"
ESSENTIAL_26G_REF="/fsimb/groups/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/proteincoding_pseudogenes_lincRNA.pipeline.gtf"
ESSENTIAL_SENSOR_REF="/fsimb/groups/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Sensor/21U_sensor_to_plot.bed"
ESSENTIAL_REPEAT_MASKER="/fsimb/groups/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/repeatMasker/c_elegans.Wb235.ce11.repeatMasker.noSimpleLow.bed"
ESSENTIAL_REPEAT_REF="/fsimb/groups/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/repEnrichNoSimpleLow"
ESSENTIAL_BIOTYPES_TABLE="/fsimb/groups/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/wormbaseGeneID2biotype.txt"

ESSENTIAL_SPECIES="worm"   // necessary for miRDeep2, used to refer to UCSC
ESSENTIAL_SAMPLE_PREFIX="Sample_"
ESSENTIAL_FEATURETYPE="gene_type" //gencode uses gene_type; ensemble uses gene_biotype
ESSENTIAL_STRANDED="reverse"
ESSENTIAL_PAIRED="no"
ESSENTIAL_THREADS=4
ESSENTIAL_MISMATCHES=0

NORMALIZATION_TO_NONSTRUCT="no"
ESSENTIAL_ORG=""           // UCSC organism
ESSENTIAL_DB="" 
REPENRICH_BED="TRUE" // the annotation was prepared as bed file.

// read size to keep (added 8 bp to account for barcodes)
ESSENTIAL_ADAPTER_SEQUENCE="TGGAATTCTCGGGTGCCAAGG"
ESSENTIAL_MINREADLENGTH=26
ESSENTIAL_MAXREADLENGTH=38
ESSENTIAL_MINADAPTEROVERLAP=5
*/


// ZEBRAFISH
ESSENTIAL_PROJECT="/project/"
ESSENTIAL_STAR_REF="/fsimb/groups/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv10/index/star/2.5.2b"
ESSENTIAL_GENESGTF="/fsimb/groups/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv10/Annotation/Genes/Danio_rerio.GRCz10.80.chr.gtf"
ESSENTIAL_RRNA_BOWTIE_REF="/fsimb/groups/imb-kettinggr/genomes/Danio_rerio/NCBI/rRNA/bowtie/1.2.2/zebrafish_all_rRNA" // necessary for fastqscreen
ESSENTIAL_GENESBED="/fsimb/groups/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv10/Annotation/Genes/Danio_rerio.GRCz10.80.chr.bed"
ESSENTIAL_BOWTIE_REF="/fsimb/groups/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv10/index/bowtie/1.2.1.1/GRCz10"
ESSENTIAL_GENOME_REF="/fsimb/groups/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv10/Sequence/WholeGenomeFastaTopLevel/genome.fa"

// vars for piRNA analyis
ESSENTIAL_FEATURES="/fsimb/groups/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv10/Annotation/Genes/pipeline/transposons.bed"
ESSENTIAL_GENES="/fsimb/groups/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv10/Annotation/Genes/pipeline/Danio_rerio_and_repeat_masker.GRCz10.80.chr.gtf"
ESSENTIAL_REPEAT_MASKER="/home/adomingu/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv10/Sequence/repeatMasker/danRer10.noSimpleLow.fa.out"
ESSENTIAL_REPEAT_REF="/home/adomingu/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv10/Sequence/repEnrichNoSimpleLow"
REPENRICH_BED="FALSE" 

ESSENTIAL_SPECIES="zebrafish"   // necessary for miRDeep2, used to refer to UCSC
ESSENTIAL_SAMPLE_PREFIX="Sample_"
ESSENTIAL_FEATURETYPE="gene_type" //gencode uses gene_type; ensemble uses gene_biotype
ESSENTIAL_DB="danRer10"
ESSENTIAL_STRANDED="yes"
ESSENTIAL_PAIRED="no"
ESSENTIAL_MISMATCHES=2
ESSENTIAL_THREADS=4

ESSENTIAL_ADAPTER_SEQUENCE="TGGAATTCTCGGGTGCCAAGG"
ESSENTIAL_READLENGTH=51      // actual read length in original raw data (incl. insert, UMIs, adapter)
ESSENTIAL_MINREADLENGTH=28
ESSENTIAL_MAXREADLENGTH=45
ESSENTIAL_MINADAPTEROVERLAP=5

// vars for mirDeep2
ESSENTIAL_MATURE_MIRNA="/fsimb/groups/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv10/Annotation/SmallRNA/mature.noSpaces.fa"
ESSENTIAL_HAIRPIN_MIRNA="~/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv10/Annotation/SmallRNA/hairpin.noSpaces.fa"


//global vars that will be reused in some global vars
PROJECT=ESSENTIAL_PROJECT
RESULTS=PROJECT + "/results"
LOGS=PROJECT + "/logs"
QC=PROJECT + "/qc"
REPORTS=PROJECT + "/reports"
PROCESSED=PROJECT + "/rawdata_processed"
MAPPED=PROJECT + "/mapped"
UNIQUE_MAPPED=MAPPED + "/unique"
UNMAPPED=MAPPED + "/unmapped"
FQ=MAPPED + "/bam2fq"
MULTI_MAPPED=MAPPED + "/multimapped"
TMP=PROJECT + "/tmp"
TRACKS=PROJECT + "/tracks"
