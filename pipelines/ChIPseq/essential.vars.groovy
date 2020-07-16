ESSENTIAL_PROJECT="/your/project/folder/"   // full project directory path
ESSENTIAL_BOWTIE_REF="/fsimb/common/genomes/saccharomyces_cerevisiae/ucsc/saccer3/canonical/index/bowtie2/genome"  // full path to reference index files (base name) needed for Bowtie2
ESSENTIAL_BOWTIE_GENOME="/fsimb/common/genomes/saccharomyces_cerevisiae/ucsc/saccer3/canonical/genome/genome.fa" // full path to the reference genome FASTA file
ESSENTIAL_CHROMSIZES="/fsimb/common/genomes/saccharomyces_cerevisiae/ucsc/saccer3/canonical/genome/genome.fa.fai" // chromosome sizes file of the reference genome
ESSENTIAL_SAMPLE_PREFIX=""  // sample prefix to be trimmed in the results and reports
ESSENTIAL_BLACKLIST=""
ESSENTIAL_BSGENOME="BSgenome.Scerevisiae.UCSC.sacCer3"  // Bioconductor genome reference used by some modules
ESSENTIAL_TXDB="TxDb.Scerevisiae.UCSC.sacCer3.sgdGene" // needed for peak annotation
ESSENTIAL_ANNODB="org.Sc.sgd.db"                    // needed for peak annotation
ESSENTIAL_READLEN=51  // read length
ESSENTIAL_FRAGLEN=200   // mean length of library inserts and also minimum peak size called by MACS2
ESSENTIAL_DUP="auto"  // how MACS2 deals with duplicated reads or fragments ("auto", "no" or 1)
ESSENTIAL_MACS2_GSIZE="10000000"  // mappable genome size for MACS2 (use approx. number or "hs" for human, "mm" for mouse, "ce" for worm, "dm" for fly)
ESSENTIAL_THREADS=4   // number of threads for parallel tasks
ESSENTIAL_DB="sacCer3"  // UCSC assembly version for GREAT analysis (only for UCSC hg19, hg18, mm9 and danRer7)
ESSENTIAL_FRAGMENT_USAGE="no" // "no" for SR data; "yes" for PE data to make bigWig tracks with reconstituted fragments
ESSENTIAL_PAIRED="no"   // to perform MACS2 peak calling in SR mode ("no") or PE mode ("yes") 
ESSENTIAL_STRANDED="no"  // library prep protocol strandness: no|yes|reverse
ESSENTIAL_BAMCOVERAGE="--binSize 10 --normalizeUsing CPM"  // deepTools options for making normalised bigWig tracks
ESSENTIAL_USE_BOWTIE1=false // bowtie1 may be used for se design with very short readlength (<50nt). By default (false) bowtie2 is used.

// project folders
PROJECT=ESSENTIAL_PROJECT
LOGS=PROJECT + "/logs"
MAPPED=PROJECT + "/mapped"
QC=PROJECT + "/qc"
REPORTS=PROJECT + "/reports"
RESULTS=PROJECT + "/results"
TMP=PROJECT + "/tmp"
TRACKS=PROJECT + "/tracks"

// optional pipeline stages to include
RUN_IN_PAIRED_END_MODE=(ESSENTIAL_PAIRED == "yes")
RUN_PEAK_ANNOTATION=true
RUN_DIFFBIND=true
RUN_TRACKHUB=false
// RUN_USING_UNFILTERED_BAM=false  // deprecated 
// Pipeline performes 2 branches for bam file processing in parallel
// default branch: remove duplicates and multimappers (MapQ filtered) from BAM file
// assumes that duplicated library fragments are unwanted PCR artifacts
// discards multimapping reads as habitually done in most ChIP-seq studies.
// unfilt branch: alternative workflow using the unfiltered BAM files
// may be preferable when studying some types of repetitive regions.
// Make sure to have ESSENTIAL_DUP="auto" for MACS2 peak calling.
