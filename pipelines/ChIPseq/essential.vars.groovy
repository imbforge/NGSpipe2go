// ChIP-Seq ESSENTIAL VARIABLES

// define essential variables here
// further module-specific variables can be adjusted in the corresponding .header files for each module
//
// The pipeline performes 2 branches for bam file processing in parallel.
// filtered branch: remove multimapping reads (MapQ filtered) as habitually done in most ChIP-seq studies
//     and (optionally) removes duplicates from BAM file (assuming that duplicated library fragments are 
//     unwanted PCR artifacts. Duplicate removal is not recommended for single end design!).
// unfiltered branch: alternative workflow using the unfiltered BAM files. This may be preferable when 
//     studying some types of repetitive regions (Make sure to have ESSENTIAL_DUP="auto" for MACS2 peak calling).

ESSENTIAL_PROJECT="/your/project/folder/"   // full project directory path
ESSENTIAL_BOWTIE_REF="/fsimb/common/genomes/saccharomyces_cerevisiae/ucsc/saccer3/canonical/index/bowtie2/genome"  // full path to reference index files (base name) needed for Bowtie2
ESSENTIAL_BOWTIE_GENOME="/fsimb/common/genomes/saccharomyces_cerevisiae/ucsc/saccer3/canonical/genome/genome.fa"  // full path to the reference genome FASTA file
ESSENTIAL_CHROMSIZES="/fsimb/common/genomes/saccharomyces_cerevisiae/ucsc/saccer3/canonical/genome/genome.fa.fai" // chromosome sizes file of the reference genome
ESSENTIAL_SAMPLE_PREFIX=""     // sample prefix to be trimmed in the results and reports
ESSENTIAL_BSGENOME="BSgenome.Scerevisiae.UCSC.sacCer3"  // Bioconductor genome reference used by some modules
ESSENTIAL_TXDB="TxDb.Scerevisiae.UCSC.sacCer3.sgdGene" // needed for peak annotation
ESSENTIAL_ANNODB="org.Sc.sgd.db"                    // needed for peak annotation
ESSENTIAL_READLEN=51           // read length
ESSENTIAL_FRAGLEN=200          // mean length of library inserts (default 200)
ESSENTIAL_MACS2_BROAD=false    // use broad setting for broad peak calling in MACS2 (default false)
ESSENTIAL_PEAK_MINLENGTH=200   // MACS2 minimum peak length. MACS2 default is fragment size. Should be increased if broad option is used.
ESSENTIAL_DEDUPLICATION=false  // Remove duplicated reads in filtered branch (default false is strongly recommended for single end data!).  
ESSENTIAL_DUP="auto"           // how MACS2 deals with duplicated reads or fragments ("auto", "all" or 1)
ESSENTIAL_MACS2_GSIZE="10000000"  // mapable genome size for MACS2 (use approx. number or "hs" for human, "mm" for mouse, "ce" for worm, "dm" for fly)
ESSENTIAL_BLACKLIST=""         // path to bed file with blacklisted regions (default: empty string). Peaks in these regions will be removed from the peakset. 
ESSENTIAL_THREADS=4            // number of threads for parallel tasks
ESSENTIAL_DB="sacCer3"         // UCSC assembly version for GREAT analysis (only for UCSC hg19, hg38, mm9 and mm10)
ESSENTIAL_FRAGMENT_USAGE="no"  // "no" for SR data; "yes" for PE data to make bigWig tracks with reconstituted fragments
ESSENTIAL_PAIRED="no"          // to perform MACS2 peak calling in SR mode ("no") or PE mode ("yes") 
ESSENTIAL_STRANDED="no"        // library prep protocol strandness: no|yes|reverse
ESSENTIAL_BAMCOVERAGE="--binSize 10 --normalizeUsing CPM"  // deepTools options for making normalised bigWig tracks
ESSENTIAL_USE_BOWTIE1=false    // bowtie1 may be used for se design with very short readlength (<50nt). By default (false) bowtie2 is used.
ESSENTIAL_FULLLIBRARYSIZE=true // use total number of reads in bam for normalization in DiffBind analysis (default true; false: use only reads overlapping peaks)

//Adapter trimming with Cutadapt (optional). Not necessary for paired end libraries.
RUN_CUTADAPT=false
ESSENTIAL_ADAPTER_SEQUENCE="Illumina=CTGTCTCTTATACACATCT" // standard sequence to trim illumina reads 
ESSENTIAL_MINADAPTEROVERLAP=3  // minimal overlap of the read and the adapter for an adapter to be found (default 3)
ESSENTIAL_MINREADLENGTH=20     // minimal length of reads to be kept


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

// optional pipeline stages to include
RUN_IN_PAIRED_END_MODE=(ESSENTIAL_PAIRED == "yes")
RUN_PEAK_ANNOTATION=true
RUN_DIFFBIND=true
RUN_ENRICHMENT=true
RUN_TRACKHUB=false


