ESSENTIAL_PROJECT="/fsimb/groups/imb-bioinfocf/projects/cfb_internal/tmp/ngspipe2go_rnaseq_test"
ESSENTIAL_STAR_REF="/fsimb/common/genomes/saccharomyces_cerevisiae/ensembl/R64/canonical/index/star/2.5.1b/"
ESSENTIAL_GENESGTF="/fsimb/common/genomes/saccharomyces_cerevisiae/ensembl/R64/canonical/annotation/Saccharomyces_cerevisiae.R64-1-1.86.gtf"
ESSENTIAL_GENESBED="/fsimb/common/genomes/saccharomyces_cerevisiae/ensembl/R64/canonical/annotation/Saccharomyces_cerevisiae.R64-1-1.86.bed"
ESSENTIAL_CHROMSIZES="/fsimb/common/genomes/saccharomyces_cerevisiae/ensembl/R64/canonical/index/star/2.5.1b/chrNameLength.txt"  // chromosome sizes file
ESSENTIAL_FEATURETYPE="gene_biotype" //gencode uses gene_type; ensemble uses gene_biotype
ESSENTIAL_SAMPLE_PREFIX=""
ESSENTIAL_PAIRED="no"           // paired end design
ESSENTIAL_STRANDED="reverse"    // strandness: no|yes|reverse
ESSENTIAL_ORG="yeast"           // UCSC organism
ESSENTIAL_DB="sacCer3"          // UCSC assembly version
ESSENTIAL_READLENGTH=50         // added for STAR version > 2.4.1a
ESSENTIAL_THREADS=4             // number of threads for parallel tasks
ESSENTIAL_FRAGMENT_USAGE="no"   //should fragments be reconstituted? should always be no for rnaseq
ESSENTIAL_FILTER_CHR=""         //chromosomes to include in post-mapping analysis.
ESSENTIAL_BAMCOVERAGE="--binSize 1 --skipNonCoveredRegions --normalizeUsing CPM" // NO smoothing should be done for RNAseq
// size of bins (in bases)  and method to normalize the number of reads per bin to generate bigwig file.

//global vars
PROJECT=ESSENTIAL_PROJECT
LOGS=PROJECT + "/logs"
MAPPED=PROJECT + "/mapped"
QC=PROJECT + "/qc"
REPORTS=PROJECT + "/reports"
RESULTS=PROJECT + "/results"
TMP=PROJECT + "/tmp"
TRACKS=PROJECT + "/tracks"
FUSION=PROJECT + "/fusion"

// optional pipeline stages to include
RUN_TRACKHUB=false
RUN_IN_PAIRED_END_MODE=(ESSENTIAL_PAIRED == "yes")
