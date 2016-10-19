ESSENTIAL_PROJECT="/project/"
ESSENTIAL_STAR_REF="/data/star_genomes/Mus_musculus/UCSC/mm9/star2.4.2a_noGTF/"
ESSENTIAL_GENESGTF="/data/igenomes_reference/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf"
ESSENTIAL_GENESBED="/data/annotation/mm9/mm9_UCSC_knownGene.bed"
ESSENTIAL_FEATURETYPE="gene_type" //gencode uses gene_type; ensemble uses gene_biotype
ESSENTIAL_SAMPLE_PREFIX="Sample_"
ESSENTIAL_PAIRED="no"           // paired end design
ESSENTIAL_STRANDED="reverse"    // strandness: no|yes|reverse
ESSENTIAL_ORG="mouse"           // UCSC organism
ESSENTIAL_DB="mm9"              // UCSC assembly version
ESSENTIAL_READLENGTH=50         // added for STAR version > 2.4.1a
ESSENTIAL_THREADS=4             // number of threads for parallel tasks

//global vars
PROJECT=ESSENTIAL_PROJECT
LOGS=PROJECT + "/logs"
MAPPED=PROJECT + "/mapped"
QC=PROJECT + "/qc"
REPORTS=PROJECT + "/reports"
RESULTS=PROJECT + "/results"
TMP=PROJECT + "/tmp"
TRACKS=PROJECT + "/tracks"
