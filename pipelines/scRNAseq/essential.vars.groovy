// Essential Variables

// General
ESSENTIAL_PROJECT="/your/project/folder/"
ESSENTIAL_SAMPLE_PREFIX="imb_gcf_2020_03_hco_20_"
ESSENTIAL_PAIRED="no"           // paired end design ("no" for MARS-Seq and 10X, because R2 in MARS_Seq and R1 in 10X contain UMI and barcodes only)
ESSENTIAL_STRANDED="yes"         // strandness: no|yes|reverse // defaults per pipeline: "yes" for tenX
ESSENTIAL_READLENGTH=130        // added for STAR version > 2.4.1a
ESSENTIAL_THREADS=4             // number of threads for parallel tasks
ESSENTIAL_FRAGMENT_USAGE="no"   //should fragments be reconstituted? should always be no for rnaseq
ESSENTIAL_FILTER_CHR=""         //chromosomes to include in post-mapping analysis.
ESSENTIAL_BAMCOVERAGE="--binSize 1 --skipNonCoveredRegions --normalizeUsing CPM" // NO smoothing should be done for RNAseq

// Alignment with STAR (used for Smart-Seq2 and MARS-Seq pipeline)
ESSENTIAL_STAR_REF="/fsimb/common/genomes/homo_sapiens/gencode/release-25_GRCh38.p7/full/index/star/2.7.3a/"
// Alignment with Cellranger (used for 10X pipeline)
ESSENTIAL_TENX_TRANSCRIPTOME="/fsimb/common/genomes/homo_sapiens/10X/grch38_ensembl98/refdata-gex-GRCh38-2020-A"
ESSENTIAL_TENX_EXPECTED_CELLS=7000 // cellranger default 
ESSENTIAL_TENX_AGGRCSV="aggregation.csv"
ESSENTIAL_TENX_NORMALIZED="mapped"

// Reference genomes for FastqScreen
ESSENTIAL_FASTQSCREEN_GENOME="Human::/fsimb/common/genomes/homo_sapiens/gencode/release-25_GRCh38.p7/full/index/bowtie2/GRCh38.p7.genome"  //bowtie2 reference for the genome the samples are from, this is used for the fastqscreen
ESSENTIAL_FASTQSCREEN=ESSENTIAL_FASTQSCREEN_GENOME + ",PHIX::/fsimb/common/genomes/phix/19930428/NCBI/index/bowtie2/2.3.4.3/ncbi_phix,ERCC::/fsimb/common/genomes/ERCC/index/bowtie2/2.3.4.3/ERCC92,rRNA::/fsimb/common/genomes/contaminants/fastqscreen_references/rrna/v1/index/bowtie2/2.3.4.3/hs_mm_ce_dm_rn_dr_xt_rRNA,Mycoplasma::/fsimb/common/genomes/contaminants/fastqscreen_references/mycoplasma/v1/index/bowtie2/2.3.4.3/mycoplasma_all_ref,E.coli::/fsimb/common/genomes/Escherichia_coli/ensembl/full/index/bowtie2/Escherichia_coli_str_k_12_substr_dh10b.ASM1942v1.31.dna.genome,B.taurus::/fsimb/common/genomes/bos_taurus/ensembl/3.1/full/index/bowtie2/2.2.9/UMD3.1" //references for fastqscreen to use if run, this are our standard references please include yours 

// Genome annotation
ESSENTIAL_GENESGTF="/fsimb/common/genomes/homo_sapiens/10X/grch38_ensembl98/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
ESSENTIAL_GENESBED="/fsimb/common/genomes/homo_sapiens/10X/grch38_ensembl98/refdata-gex-GRCh38-2020-A/genes/genes.bed"
ESSENTIAL_CHROMSIZES="/fsimb/common/genomes/homo_sapiens/10X/grch38_ensembl98/refdata-gex-GRCh38-2020-A/star/chrNameLength.txt"  // chromosome sizes file
ESSENTIAL_FEATURETYPE="gene_type" //gencode (also 10X Genomics) uses gene_type; ensemble uses gene_biotype
ESSENTIAL_ORG="human"           // UCSC organism
ESSENTIAL_DB="hg38"              // UCSC assembly version
ESSENTIAL_MTGENES=""  // list with gene_id of mitochondrial genes (give path within ESSENTIAL_PROJECT). If not given (empty string), all genes annotated on the mitochondrial chromosome are used by default.


// Barcodes (used for UMI and barcode extraction in MARS-Seq)
ESSENTIAL_WHITELIST="single_cell_barcodes_whitelist.csv" //the barcode list should be a list of valid barcodes separated by newline
ESSENTIAL_BCPATTERN="CCCCCCCNNNNNNNN" //barcode pattern as it is present in MARS-Seq data.
                                     //C stands for Cellbarcode and N is part of the random barcode
                                     //Please be aware that the addumibarcode module which uses this variable will only 
                                     //work for MARS-Seq where the barcode and the umi are present in the second reads otherwise it has to be modified! 

//Adapter trimming
RUN_CUTADAPT=false
//ESSENTIAL_ADAPTER_SEQUENCE="TruSeqLTHT=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" 
ESSENTIAL_ADAPTER_SEQUENCE="Nextera=CTGTCTCTTATACACATCT" //standard sequence to trim illumina reads
ESSENTIAL_MINADAPTEROVERLAP=5
ESSENTIAL_MINREADLENGTH=30

//global vars
PROJECT=ESSENTIAL_PROJECT
LOGS=PROJECT + "/logs"
MAPPED=PROJECT + "/mapped"
TRIMMED=PROJECT + "/trimmed"
QC=PROJECT + "/qc"
REPORTS=PROJECT + "/reports"
RESULTS=PROJECT + "/results"
TMP=PROJECT + "/tmp"
TRACKS=PROJECT + "/tracks"
FUSION=PROJECT + "/fusion"

// optional pipeline stages to include
RUN_TRACKHUB=false
RUN_IN_PAIRED_END_MODE=(ESSENTIAL_PAIRED == "yes")


