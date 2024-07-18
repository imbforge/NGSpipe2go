// ChIP-Seq ESSENTIAL VARIABLES

// Define essential variables here.
// Further module-specific variables can be adjusted in the corresponding ".header" files for each module.
//
// The pipeline performes 2 branches for BAM file processing in parallel.
// Filtered branch: removes multimapping reads (MapQ filtered) as habitually done in most ChIP-seq studies
//     and (optionally) removes duplicates from BAM file assuming they are unwanted PCR artifacts (duplicate 
//     removal is usually recommended for paired-end data but not for single-end data).
// Unfiltered branch: alternative workflow using the unfiltered BAM files. This may be preferable when 
//     studying some types of repetitive regions (keep ESSENTIAL_DUP="auto" for MACS2 peak calling).

// General parameters
ESSENTIAL_PROJECT="project_dir"   // full project directory path
ESSENTIAL_THREADS=4            // number of threads for parallel tasks
ESSENTIAL_SAMPLE_PREFIX=""     // sample name prefix to be trimmed in the results and reports

// Mapping and filtering parameters
ESSENTIAL_BOWTIE_REF="index_path/genome"  // full path to reference index files (base name) needed for Bowtie2
ESSENTIAL_BOWTIE_GENOME="genome_path/genome.fa"  // full path to the reference genome FASTA file
ESSENTIAL_CHROMSIZES="genome_path/genome.fa.fai" // chromosome sizes file of the reference genome
ESSENTIAL_PAIRED="yes"          // perform MACS2 peak calling in single-end mode ("no") or paired-end mode ("yes") according to the sequencing mode.
ESSENTIAL_STRANDED="no"        // library prep protocol strandness: no|yes|reverse (typically "no" in ChIP-seq)
ESSENTIAL_READLEN=42           // read length in the FastQ files
ESSENTIAL_FRAGLEN=200          // mean length of library inserts (typically around 200 bp in ChIP-seq but may be lower in DIP-seq, CUT&RUN etc. libraries).
ESSENTIAL_FRAGMENT_USAGE="yes" // should fragments be reconstituted for generating bigWig coverage tracks? single reads are extended to ESSENTIAL_FRAGLEN, paired reads are extended to match the fragment size defined by the read mates.
ESSENTIAL_DEDUPLICATION=(ESSENTIAL_PAIRED == "yes") // remove duplicated reads in the filtered branch of the pipeline for paired-end but not for single-end data (may need to be changed in case of ultra-deep PE sequencing of small genomes).  
ESSENTIAL_BAMCOVERAGE="--binSize 10 --normalizeUsing CPM"  // deepTools options for making normalised bigWig tracks

// Annotation parameters
ESSENTIAL_BSGENOME="BSgenome.Scerevisiae.UCSC.sacCer3"  // Bioconductor genome reference used by some modules
ESSENTIAL_TXDB="TxDb.Scerevisiae.UCSC.sacCer3.sgdGene"  // needed for peak annotation
ESSENTIAL_ANNODB="org.Sc.sgd.db"  // needed for peak annotation
ESSENTIAL_DB="sacCer3"            // UCSC assembly version for GREAT analysis (only for UCSC hg19, hg38, mm9 and mm10)
ESSENTIAL_BLACKLIST=""           // path to a BED file with blacklisted regions (default: empty string). Peaks in these regions will be removed from the peakset. IMPORTANT NOTE: blacklist is skipped if RUN_MAKE_GREYLIST is set to true below! 

// FASTQ-Screen parameters
ESSENTIAL_FASTQSCREEN_PERC=1    // contaminant filter, if a contaminant is consuming at least this percentage of reads in at least one sample, contaminant will be shown in report
ESSENTIAL_FASTQSCREEN_GENOME="Yeast::/fsimb/common/genomes/saccharomyces_cerevisiae/ensembl/R64/canonical/index/bowtie2/2.3.2/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel"  //bowtie2 reference for the genome the samples are from, this is used for the fastqscreen
ESSENTIAL_FASTQSCREEN=ESSENTIAL_FASTQSCREEN_GENOME + ",PHIX::/fsimb/common/genomes/phix/19930428/NCBI/index/bowtie2/2.3.4.3/ncbi_phix,ERCC::/fsimb/common/genomes/ERCC/index/bowtie2/2.3.4.3/ERCC92,rRNA::/fsimb/common/genomes/contaminants/fastqscreen_references/rrna/v1/index/bowtie2/2.3.4.3/hs_mm_ce_dm_rn_dr_xt_rRNA,Mycoplasma::/fsimb/common/genomes/contaminants/fastqscreen_references/mycoplasma/v1/index/bowtie2/2.3.4.3/mycoplasma_all_ref,E.coli::/fsimb/common/genomes/Escherichia_coli/ensembl/full/index/bowtie2/Escherichia_coli_str_k_12_substr_dh10b.ASM1942v1.31.dna.genome,B.taurus::/fsimb/common/genomes/bos_taurus/ensembl/3.1/full/index/bowtie2/2.2.9/UMD3.1" //references for fastqscreen to use if run, this are our standard references please include yours 


// Adapter trimming with Cutadapt (optional). Usually not needed when the reads are much shorter than the library inserts.
// Cutadapt recommends using full length adapter sequences since adapter fragments might occur in the genome
RUN_CUTADAPT=false
ESSENTIAL_ADAPTER_SEQUENCE="Illumina=AGATCGGAAGAG"
ESSENTIAL_ADAPTER_SEQUENCE_R2="" // in case of paired end sequencing the R2 adapter sequence (which will be trimmed from the 3' end of R2, -A argument in Cutadapt), by default the same sequence as defined in ESSENTIAL_ADAPTER_SEQUENCE is used
ESSENTIAL_MINADAPTEROVERLAP=3  // minimal overlap of the read and the adapter for an adapter to be found (default 3)
ESSENTIAL_MINREADLENGTH=20     // minimal length of reads to be kept
ESSENTIAL_BASEQUALCUTOFF=20    // trim low-quality ends from reads (if nextseqtrim is true, qualities of terminal G bases are ignored)  
ESSENTIAL_NEXTSEQTRIM=true     // accounts for terminal G bases during base quality trimming incorporated by faulty dark cycles observed with two-color chemistry (as in NextSeq) 

// Peak calling with MACS2
ESSENTIAL_MACS2_BROAD=false    // use "true" for broad peak calling in MACS2 (default: "false")
ESSENTIAL_DUP="auto"           // how MACS2 deals with duplicated reads or fragments: "auto" (default), "all" or 1
ESSENTIAL_MACS2_GSIZE="10000000"  // mapable genome size for MACS2 (approx. size in bp or use "hs" for human, "mm" for mouse, "ce" for worm, "dm" for fly)
ESSENTIAL_MIN_PEAKLENGTH=""  // MACS2 minimum peak length. Default (empty string) uses predicted fragment size. Could be increased if broad option is used. For ATAC-Seq reduce to 100.

// Differential binding analysis with DiffBind
// Note that DiffBind works with "default" parameters which depend on the 
// context of other parameter settings (see DiffBind documentation for 
// explanation). Since version 3, DiffBind includes the data from ALL samples in 
// a single model. 
// In this pipeline, you can specify sub_experiments in contrasts_diffbind.txt
// to define those sample groups which shall be combined in a model.
RUN_DIFFBIND=true
ESSENTIAL_DIFFBIND_LIBRARY="default" // DiffBind method to calculate library sizes. One of "full", "RiP", "background" and "default"  
ESSENTIAL_DIFFBIND_NORM="default"    // DiffBind method to calculate normalization factors. One of "lib", "RLE", "TMM", "native" and "default". 
ESSENTIAL_SUMMITS=200                // Re-center peaks around consensus summit with peak width 2x ESSENTIAL_SUMMITS (0 means no re-centering)

// further optional pipeline stages to include
RUN_IN_PAIRED_END_MODE=(ESSENTIAL_PAIRED == "yes")
RUN_FASTQSCREEN=true            // check for contaminations using FastQ Screen
RUN_MAKE_GREYLIST=true          // greylist is generated from all Control files and is then applied to MACS2 peak files like a blacklist (default: true). IMPORTANT NOTE: If true, a bed file given in ESSENTIAL_BLACKLIST is ignored! 
RUN_PEAK_ANNOTATION=true
RUN_ENRICHMENT=true
RUN_TRACKHUB=false
RUN_UPSETPLOT=false // For larger projects, the calculation of the respective combination matrix may take several hours.

// project folders
PROJECT=ESSENTIAL_PROJECT
LOGS=PROJECT + "/logs"
MAPPED=PROJECT + "/mapped"
TRIMMED=PROJECT + "/trimmed"
QC=PROJECT + "/qc"
REPORTS=PROJECT + "/reports"
RESULTS=PROJECT + "/results"
TMP=PROJECT + "/tmp"
TRACKS=PROJECT + "/tracks"

