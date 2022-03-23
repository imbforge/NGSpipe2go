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
ESSENTIAL_BLACKLIST=""            // path to a BED file with blacklisted regions (default: empty string). Peaks in these regions will be removed from the peakset. 

// Adapter trimming with Cutadapt (optional). Usually not needed when the reads are much shorter than the library inserts.
RUN_CUTADAPT=false
ESSENTIAL_ADAPTER_SEQUENCE="Illumina=CTGTCTCTTATACACATCT" // standard sequence to trim illumina reads 
ESSENTIAL_MINADAPTEROVERLAP=3  // minimal overlap of the read and the adapter for an adapter to be found (default 3)
ESSENTIAL_MINREADLENGTH=20     // minimal length of reads to be kept
ESSENTIAL_BASEQUALCUTOFF=20    // trim low-quality ends from reads (if nextseqtrim is true, qualities of terminal G bases are ignored)  
ESSENTIAL_NEXTSEQTRIM=true     // accounts for terminal G bases during base quality trimming incorporated by faulty dark cycles observed with two-color chemistry (as in NextSeq) 

// Peak calling with MACS2
ESSENTIAL_MACS2_BROAD=false    // use "true" for broad peak calling in MACS2 (default: "false")
ESSENTIAL_DUP="auto"           // how MACS2 deals with duplicated reads or fragments: "auto" (default), "all" or 1
ESSENTIAL_MACS2_GSIZE="10000000"  // mapable genome size for MACS2 (approx. size in bp or use "hs" for human, "mm" for mouse, "ce" for worm, "dm" for fly)

// Differential binding analysis with DiffBind
// Note that DiffBind3 works with "default" parameters which depend on the 
// context of other parameter settings (see DiffBind documentation for 
// explanation). Unlike DiffBind2, DiffBind3 includes the data from ALL samples in 
// a single model. 
// In this pipeline, you can specify sub_experiments in contrasts_diffbind.txt
// to define those sample groups which shall be combined in a model.
RUN_DIFFBIND=true
ESSENTIAL_DIFFBIND_VERSION=3         // Beginning with version 3, DiffBind has included new functionalities and modified default settings. Earlier versions are also supported here.
ESSENTIAL_DIFFBIND_LIBRARY="default" // DiffBind method to calculate library sizes. One of "full", "RiP", "background" and "default"  
ESSENTIAL_DIFFBIND_NORM="default"    // DiffBind method to calculate normalization factors. One of "lib", "RLE", "TMM", "native" and "default". Not applicable for DiffBind2.

// further optional pipeline stages to include
RUN_IN_PAIRED_END_MODE=(ESSENTIAL_PAIRED == "yes")
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

