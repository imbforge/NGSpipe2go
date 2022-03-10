// RNA-Seq ESSENTIAL VARIABLES

// Define essential variables here.
// Further module-specific variables can be adjusted in the corresponding ".header" files for each module.

// General parameters
ESSENTIAL_PROJECT="/project/"   // full project directory path
ESSENTIAL_THREADS=4             // number of threads for parallel tasks
ESSENTIAL_SAMPLE_PREFIX=""      // sample name prefix to be trimmed in the results and reports

// Mapping parameters
ESSENTIAL_BOWTIE_REF="index_path/genome"  // full path to reference index files (base name) needed for Bowtie2
ESSENTIAL_BOWTIE_GENOME="genome_path/genome.fa"  // full path to the reference genome FASTA file
ESSENTIAL_CHROMSIZES="genome_path/genome.fa.fai" // chromosome sizes file of the reference genome
ESSENTIAL_PAIRED="yes"          // single-end ("no") or paired-end ("yes") data
ESSENTIAL_READLEN=50            // read length in the FastQ files
ESSENTIAL_FRAGLEN=500           // mean length of library inserts
ESSENTIAL_STRANDED="no"         // library prep protocol strandness: no|yes|reverse
ESSENTIAL_FILTER_CHR=""         // chromosomes to include in post-mapping analysis
ESSENTIAL_BAMCOVERAGE="--binSize 1 --skipNonCoveredRegions --normalizeUsing CPM"  // size of bins (in bp) and normalisation method for the bigWig tracks

// Annotation parameters
ESSENTIAL_ORG="mouse"           // UCSC organism name (human, mouse, fly, worm, yeast, zebrafish)
ESSENTIAL_DB="mm10"             // UCSC assembly version (GREAT analysis only for UCSC hg19, hg38, mm9 and mm10)
ESSENTIAL_BSGENOME="BSgenome.Mmusculus.UCSC.mm10"   // Bioconductor genome reference used by some modules
ESSENTIAL_TXDB="TxDb.Mmusculus.UCSC.mm10.knownGene"  // needed for peak annotation
ESSENTIAL_ANNODB="org.Mm.eg.db"   // needed for peak annotation
ESSENTIAL_BLACKLIST=""            // path to a BED file with blacklisted regions (default: empty string). Peaks in these regions will be removed from the peakset.


// Capture regions files (only relevant for CapSTARR-seq; use empty strings for conventional STARR-seq)
ESSENTIAL_GENESBED=ESSENTIAL_PROJECT+"/capture_regions.bed"  // for CapSTARR-seq: capture regions BED file
ESSENTIAL_GENESGTF=ESSENTIAL_PROJECT+"/capture_regions.gtf"  // for CapSTARR-seq: capture regions GTF file (each capture region should be a single-exon, single transcript gene with "gene_id" and "gene_name" fields. Only "exon" GTF entries are required)

// Peak calling with STARRPeaker (can be used for both conventional STARR-seq and CapSTARR-seq)
RUN_STARRPEAKER=false
ESSENTIAL_STARRPEAKER_BINBED="/full/path/to/genomic_bins.bin.bed"  // BED file with genomic bins used in peak calling (should end in ".bin.bed")
ESSENTIAL_STARRPEAKER_COVTSV="/full/path/to/genomic_bin_covariates.tsv"   // tab-separated file with covariate scores per genomic bin (one covariate per column). This file can be created from the two variables below with the STARRPeaker_preproc.pipeline.groovy pipeline
ESSENTIAL_STARRPEAKER_COVBIGWIGS="/path/to/covariate1.bw,/path/to/covariate2.bw"   // comma-separated list if bigwig files to generate covariates file from
ESSENTIAL_STARRPEAKER_LINFOLDTSV="/full/path/to/linearfold_output.tsv"   // LinearFold output file
ESSENTIAL_STARRPEAKER_LINFOLDEXE="/full/path/to/linearfold_v"            // full path to LinearFold executable


// STARR-seq specific parameters
CAPSTARRSEQ_DIFFEXP=false       // for CapSTARR-seq: whether to use differential expression or fold change approach. Differential expression requires the STARR-seq DNA to have been extracted from the samples, as well as the RNA


// Optional pipeline stages to include
RUN_TRACKHUB=false              // prepare a Track Hub for the UCSC genome browser
RUN_FASTQSCREEN=true            // check for contaminations using FastQ Screen
RUN_CUTADAPT=false              // optional read trimming with Cutadapt e.g. if using high read length
RUN_IN_PAIRED_END_MODE=(ESSENTIAL_PAIRED == "yes") // no need to change this line
// Optional pipeline stages specific to conventional STARR-seq (not CapSTARR-seq)
RUN_PEAK_ANNOTATION=true
RUN_ENRICHMENT=true
RUN_UPSETPLOT=false // For larger projects, the calculation of the respective combination matrix may take several hours.


// FASTQ-Screen parameters
ESSENTIAL_FASTQSCREEN_PERC=1    // contaminant filter, if a contaminant is consuming at least this percentage of reads in at least one sample, contaminant will be shown in report
ESSENTIAL_FASTQSCREEN_GENOME="Yeast::/fsimb/common/genomes/saccharomyces_cerevisiae/ensembl/R64/canonical/index/bowtie2/2.3.2/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel"  //bowtie2 reference for the genome the samples are from, this is used for the fastqscreen
ESSENTIAL_FASTQSCREEN=ESSENTIAL_FASTQSCREEN_GENOME + ",PHIX::/fsimb/common/genomes/phix/19930428/NCBI/index/bowtie2/2.3.4.3/ncbi_phix,ERCC::/fsimb/common/genomes/ERCC/index/bowtie2/2.3.4.3/ERCC92,rRNA::/fsimb/common/genomes/contaminants/fastqscreen_references/rrna/v1/index/bowtie2/2.3.4.3/hs_mm_ce_dm_rn_dr_xt_rRNA,Mycoplasma::/fsimb/common/genomes/contaminants/fastqscreen_references/mycoplasma/v1/index/bowtie2/2.3.4.3/mycoplasma_all_ref,E.coli::/fsimb/common/genomes/Escherichia_coli/ensembl/full/index/bowtie2/Escherichia_coli_str_k_12_substr_dh10b.ASM1942v1.31.dna.genome,B.taurus::/fsimb/common/genomes/bos_taurus/ensembl/3.1/full/index/bowtie2/2.2.9/UMD3.1" //references for fastqscreen to use if run, this are our standard references please include yours 

// Adapter trimming with Cutadapt (additional adapter sequences for R1 and/or R2 can be specified in the cutadapt.header file)
ESSENTIAL_ADAPTER_SEQUENCE="Nextera=CTGTCTCTTATACACATCT" // 3'adapter sequence for R1. 
ESSENTIAL_MINADAPTEROVERLAP=5
ESSENTIAL_MINREADLENGTH=30
ESSENTIAL_NEXTSEQTRIM=true   // accounts for terminal G bases during base quality trimming incorporated by faulty dark cycles observed with two-color chemistry (as in NextSeq) 


// Peak calling with MACS2 (only required for conventional STARR-seq)
ESSENTIAL_MACS2_BROAD=false    // use "true" for broad peak calling in MACS2 (default: "false")
ESSENTIAL_DUP="all"           // how MACS2 deals with duplicated reads or fragments: "auto" (default), "all" or 1
ESSENTIAL_MACS2_GSIZE="mm"  // mapable genome size for MACS2 (approx. size in bp or use "hs" for human, "mm" for mouse, "ce" for worm, "dm" for fly)

// Differential binding analysis with DiffBind (only required for conventional STARR-seq)
// Note that DiffBind3 works with "default" parameters which depend on the 
// context of other parameter settings (see DiffBind documentation for 
// explanation). Unlike DiffBind2, DiffBind3 includes the data from ALL samples in 
// a single model.
// **IMPORTANT NOTE:** Mind that for DiffBind 3 the contrasts are processed 
// isolated from each other. This means that a consensus peakset is generated 
// separately for each contrast and does not incorporate peaks called in other 
// pulldowns unrelated to the contrast. A pipeline update that allows processing 
// together pulldowns from different contrasts is currently under construction.
RUN_DIFFBIND=true
ESSENTIAL_DIFFBIND_VERSION=3         // Beginning with version 3, DiffBind has included new functionalities and modified default settings. Earlier versions are also supported here.
ESSENTIAL_DIFFBIND_LIBRARY="default" // DiffBind method to calculate library sizes. One of "full", "RiP", "background" and "default"  
ESSENTIAL_DIFFBIND_NORM="default"    // DiffBind method to calculate normalization factors. One of "lib", "RLE", "TMM", "native" and "default". Not applicable for DiffBind2.


// DESeq2 specific parameters (only required for CapSTARR-seq)
ESSENTIAL_DESEQ2_FDR=0.01       // FDR significance cutoff in the DESeq2 model (may be increased for noisy or underpowered datasets)
ESSENTIAL_DESEQ2_FC=1           // optional Fold-Change cutoff to incorporate into the DESeq2 model (default "FC=1" means no Fold-Change filter is used)
// Using FC threshold in the DESeq2 model is usually more conservative than post-hoc gene filtering by FC (which should anyway be avoided, see https://doi.org/10.1093/bib/bbab053) 


//project folders
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
