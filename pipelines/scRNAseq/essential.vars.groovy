// Essential Variables

// General
ESSENTIAL_PROJECT="project_dir"  // full project directory path
ESSENTIAL_SEQTYPE="tenX"         // sequencing type, one of "tenX", "tenXmultiome", "ParseBio", "ScaleBio", "SmartSeq"
ESSENTIAL_SAMPLE_PREFIX=""      // sample name prefix to be trimmed in the results and reports
ESSENTIAL_THREADS=16             // number of threads for parallel tasks
ESSENTIAL_PAIRED="no"           // paired end design ("no" for MARS-Seq and 10X, because R2 in MARS_Seq and R1 in 10X contain UMI and barcodes only)
ESSENTIAL_STRANDED="yes"         // strandness: no|yes|reverse // defaults per pipeline: "yes" for tenX
ESSENTIAL_FRAGLEN=200          // mean length of library inserts (default 200)
ESSENTIAL_FRAGMENT_USAGE="no"  // "no" for SR data; "yes" for PE data to make bigWig tracks with reconstituted fragments
ESSENTIAL_BAMCOVERAGE="--binSize 1 --skipNonCoveredRegions --normalizeUsing CPM" // NO smoothing should be done for RNAseq

// path to reference genome as requested for the respective sequencing assay
ESSENTIAL_GENOME_REFERENCE="/fsimb/common/genomes/homo_sapiens/10X/grch38_ensembl98/refdata-gex-GRCh38-2020-A"
ESSENTIAL_EXPECTED_CELLS=7000 // cellranger default 
ESSENTIAL_NUCLEI="no"         // set to "yes" if run with nuclei instead of cells (e.g. the RNA data of 10X Multiome assay)
ESSENTIAL_USE_AGGR_DATA=true // if true (default) load aggregated sample data for downstream processing. Otherwise, load individual sample data and aggregate in Seurat.

// Genome annotation
ESSENTIAL_GENESGTF="/fsimb/common/genomes/homo_sapiens/10X/grch38_ensembl98/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
ESSENTIAL_GENESBED="/fsimb/common/genomes/homo_sapiens/10X/grch38_ensembl98/refdata-gex-GRCh38-2020-A/genes/genes.bed"
ESSENTIAL_CHROMSIZES="/fsimb/common/genomes/homo_sapiens/10X/grch38_ensembl98/refdata-gex-GRCh38-2020-A/star/chrNameLength.txt"  // chromosome sizes file
ESSENTIAL_FEATURETYPE="gene_type" //gencode (also 10X Genomics) uses gene_type; ensemble uses gene_biotype
ESSENTIAL_ORG="human"           // UCSC organism
ESSENTIAL_DB="hg38"              // UCSC assembly version
ESSENTIAL_MTGENES="MT-"  // prefix of mitochondrial gene ids or text file with gene ids (give path within ESSENTIAL_PROJECT)
ESSENTIAL_CELLTYPE_ANNO = ["Marker"] // select celltype annotation method (one or more of "Marker", "Seurat"). The first in the list is used for downstream processing. 
// if "Marker" selected, specify marker table in CTannoMarker.header. For "Seurat" specify reference dataset in CTannoSeurat.header.

// Optional for 10X pipelines: sample demultiplexing if sample pooling per GEM well was applied (specify one line per demultiplexed sample in targets.txt)
// Either "demux_GT" for demultiplexing by genetic variance, "demux_GT_noAssignment" for demultiplexing by genetic variance without subsequent assignment of corresponding individuals across files, 
// "demux_HTO" for cell hashing or empty string for no demultiplexing (default).
// if "demux_HTO" provide "file_HTO" (fastq file with corresponding HTO sequences), "seq_HTO" (HTO sequence) and "name_HTO" (HTO name) in targets.txt 
RUN_DEMUX=""  

// Reference genomes for FastqScreen
ESSENTIAL_FASTQSCREEN_PERC=1    // contaminant filter, if a contaminant is consuming at least this percentage of reads in at least one sample, contaminant will be shown in report
ESSENTIAL_FASTQSCREEN_GENOME="Human::/fsimb/common/genomes/homo_sapiens/gencode/release-25_GRCh38.p7/full/index/bowtie2/GRCh38.p7.genome"  // bowtie2 reference index files (base name) for the genome the samples are from
ESSENTIAL_FASTQSCREEN=ESSENTIAL_FASTQSCREEN_GENOME + ",PHIX::/fsimb/common/genomes/phix/19930428/NCBI/index/bowtie2/2.3.4.3/ncbi_phix,ERCC::/fsimb/common/genomes/ERCC/index/bowtie2/2.3.4.3/ERCC92,rRNA::/fsimb/common/genomes/contaminants/fastqscreen_references/rrna/v1/index/bowtie2/2.3.4.3/hs_mm_ce_dm_rn_dr_xt_rRNA,Mycoplasma::/fsimb/common/genomes/contaminants/fastqscreen_references/mycoplasma/v1/index/bowtie2/2.3.4.3/mycoplasma_all_ref,E.coli::/fsimb/common/genomes/Escherichia_coli/ensembl/full/index/bowtie2/Escherichia_coli_str_k_12_substr_dh10b.ASM1942v1.31.dna.genome,B.taurus::/fsimb/common/genomes/bos_taurus/ensembl/3.1/full/index/bowtie2/2.2.9/UMD3.1"  // bowtie2 reference index files (base name) for common sample contaminants, as well as for PhiX and ERCC library spike-ins

// Optional Adapter trimming (usually not necessary because done by manufacturer count software) 
// Cutadapt recommends using full length adapter sequences since adapter fragments might occur in the genome
RUN_CUTADAPT=false
ESSENTIAL_ADAPTER_SEQUENCE="Nextera=CTGTCTCTTATACACATCT" //standard sequence to trim illumina reads
ESSENTIAL_ADAPTER_SEQUENCE_R2="" // in case of paired end sequencing the R2 adapter sequence (which will be trimmed from the 3' end of R2, -A argument in Cutadapt), by default the same sequence as defined in ESSENTIAL_ADAPTER_SEQUENCE is used
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
RUN_BATCHCORRECT=false // if true, please remember to set a batch-variable in the sc_integrateRNA.header file. The batch-variable has to be a column name from targets file.

