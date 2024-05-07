// Essential Variables

// General
ESSENTIAL_PROJECT="/fsimb/groups/imb-bioinfocf/projects/cfb_internal/siva/alternate_module_aggr"
ESSENTIAL_SAMPLE_PREFIX=""
ESSENTIAL_PAIRED="no"           // paired end design ("no" for MARS-Seq and 10X, because R2 in MARS_Seq and R1 in 10X contain UMI and barcodes only)
ESSENTIAL_STRANDED="yes"         // strandness: no|yes|reverse // defaults per pipeline: "yes" for tenX
ESSENTIAL_READLENGTH=130        // added for STAR version > 2.4.1a
ESSENTIAL_THREADS=16             // number of threads for parallel tasks
ESSENTIAL_FRAGMENT_USAGE="no"   //should fragments be reconstituted? should always be no for rnaseq
ESSENTIAL_FILTER_CHR=""         //chromosomes to include in post-mapping analysis.
ESSENTIAL_FRAGLEN=200          // mean length of library inserts (default 200)
ESSENTIAL_FRAGMENT_USAGE="no"  // "no" for SR data; "yes" for PE data to make bigWig tracks with reconstituted fragments
ESSENTIAL_BAMCOVERAGE="--binSize 1 --skipNonCoveredRegions --normalizeUsing CPM" // NO smoothing should be done for RNAseq

// Reference for alignment with STAR (used for Smart-Seq2 and MARS-Seq pipeline)
ESSENTIAL_STAR_REF="/fsimb/common/genomes/homo_sapiens/gencode/release-25_GRCh38.p7/full/index/star/2.7.3a/"
// Reference for alignment with Cellranger (used for 10X GEX pipeline)
ESSENTIAL_TENX_TRANSCRIPTOME="/fsimb/common/genomes/homo_sapiens/10X/grch38_ensembl98/refdata-gex-GRCh38-2020-A"
// Reference for alignment with Cellranger ATAC & Cellranger ARC (used for 10X ATAC & Multiome pipelines)
ESSENTIAL_TENX_REFERENCE="/fsimb/common/genomes/homo_sapiens/10X/grch38_ensembl98/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"
ESSENTIAL_TENX_EXPECTED_CELLS=7000 // cellranger default 
ESSENTIAL_TENX_FASTQDIR=ESSENTIAL_PROJECT + "/rawdata"
ESSENTIAL_TENX_AGGRCSV="aggregation.csv"
ESSENTIAL_USE_AGGR_DATA = false // if true load cellranger aggr results. Otherwise, load individual sample data and aggregate in Seurat
ESSENTIAL_TENX_NORMALIZED="mapped"
ESSENTIAL_TENX_NUCLEI="yes"         // set to "yes" if 10X run with nuclei instead of cells

// Optional sample demultiplexing if sample pooling per GEM well was applied in 10X pipelines (specify one line per demultiplexed sample in targets.txt)
// Either "demux_HTO" for cell hashing, "demux_GT" for demultiplexing by genetic variance or empty string (no demultiplexing)
// if "demux_HTO" provide "file_HTO" (fastq file with corresponding HTO sequences), "seq_HTO" (HTO sequence) and "name_HTO" (HTO name) in targets.txt 
RUN_DEMUX="demux_GT"  

// Reference genomes for FastqScreen
ESSENTIAL_FASTQSCREEN_PERC=1    // contaminant filter, if a contaminant is consuming at least this percentage of reads in at least one sample, contaminant will be shown in report
ESSENTIAL_FASTQSCREEN_GENOME="Human::/fsimb/common/genomes/homo_sapiens/gencode/release-25_GRCh38.p7/full/index/bowtie2/GRCh38.p7.genome"  //bowtie2 reference for the genome the samples are from, this is used for the fastqscreen
ESSENTIAL_FASTQSCREEN=ESSENTIAL_FASTQSCREEN_GENOME + ",PHIX::/fsimb/common/genomes/phix/19930428/NCBI/index/bowtie2/2.3.4.3/ncbi_phix,ERCC::/fsimb/common/genomes/ERCC/index/bowtie2/2.3.4.3/ERCC92,rRNA::/fsimb/common/genomes/contaminants/fastqscreen_references/rrna/v1/index/bowtie2/2.3.4.3/hs_mm_ce_dm_rn_dr_xt_rRNA,Mycoplasma::/fsimb/common/genomes/contaminants/fastqscreen_references/mycoplasma/v1/index/bowtie2/2.3.4.3/mycoplasma_all_ref,E.coli::/fsimb/common/genomes/Escherichia_coli/ensembl/full/index/bowtie2/Escherichia_coli_str_k_12_substr_dh10b.ASM1942v1.31.dna.genome,B.taurus::/fsimb/common/genomes/bos_taurus/ensembl/3.1/full/index/bowtie2/2.2.9/UMD3.1" //references for fastqscreen to use if run, this are our standard references please include yours 

// Genome annotation
ESSENTIAL_GENESGTF="/fsimb/common/genomes/homo_sapiens/10X/grch38_ensembl98/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
ESSENTIAL_GENESBED="/fsimb/common/genomes/homo_sapiens/10X/grch38_ensembl98/refdata-gex-GRCh38-2020-A/genes/genes.bed"
ESSENTIAL_CHROMSIZES="/fsimb/common/genomes/homo_sapiens/10X/grch38_ensembl98/refdata-gex-GRCh38-2020-A/star/chrNameLength.txt"  // chromosome sizes file
ESSENTIAL_FEATURETYPE="gene_type" //gencode (also 10X Genomics) uses gene_type; ensemble uses gene_biotype
ESSENTIAL_ORG="human"           // UCSC organism
ESSENTIAL_DB="hg38"              // UCSC assembly version
ESSENTIAL_MTGENES="MT-"  // list with gene_id of mitochondrial genes (give path within ESSENTIAL_PROJECT)
ESSENTIAL_CELLTYPE_ANNO = ["Marker"] // select celltype annotation method (one or more of "Seurat", "Marker"). The first in the list is used for downstream processing. 
// if "Seurat" selected, specify reference dataset in CTannoSeurat.header. For "Marker" specify marker table in CTannoMarker.header.


// Barcodes (used for UMI and barcode extraction in MARS-Seq)
ESSENTIAL_WHITELIST="single_cell_barcodes_whitelist.csv" //the barcode list should be a list of valid barcodes separated by newline
ESSENTIAL_BCPATTERN="CCCCCCCNNNNNNNN" //barcode pattern as it is present in MARS-Seq data.
                                     //C stands for Cellbarcode and N is part of the random barcode
                                     //Please be aware that the addumibarcode module which uses this variable will only 
                                     //work for MARS-Seq where the barcode and the umi are present in the second reads otherwise it has to be modified! 

//Adapter trimming
// Cutadapt recommends using full length adapter sequences since adapter fragments might occur in the genome
RUN_CUTADAPT=false
//ESSENTIAL_ADAPTER_SEQUENCE="TruSeqLTHT=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" 
ESSENTIAL_ADAPTER_SEQUENCE="Nextera=CTGTCTCTTATACACATCT" //standard sequence to trim illumina reads
ESSENTIAL_MINADAPTEROVERLAP=5
ESSENTIAL_MINREADLENGTH=30
ESSENTIAL_BASEQUALCUTOFF=20    // trim low-quality ends from reads (if nextseqtrim is true, qualities of terminal G bases are ignored)  
ESSENTIAL_NEXTSEQTRIM=true     // accounts for terminal G bases during base quality trimming incorporated by faulty dark cycles observed with two-color chemistry (as in NextSeq) 

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
RUN_BATCHCORRECT=true // if true, please remember to set a batch-variable in the sc_integrateRNA.header file. The batch-variable has to be a column name from targets file.

