// Essential Variables

// General
ESSENTIAL_PROJECT="/your/project/folder/"
ESSENTIAL_SAMPLE_PREFIX=""
ESSENTIAL_PAIRED="yes"         // paired end, because 10X library
ESSENTIAL_STRANDED="yes"       // strandness: no|yes|reverse // defaults per pipeline: "yes" for tenX
ESSENTIAL_READLENGTH=130       // added for STAR version > 2.4.1a
ESSENTIAL_THREADS=4            // number of threads for parallel tasks
ESSENTIAL_FRAGMENT_USAGE="no"  //should fragments be reconstituted? should always be no for rnaseq
ESSENTIAL_FILTER_CHR=""        //chromosomes to include in post-mapping analysis.
ESSENTIAL_FRAGLEN=300          // mean length of library inserts (default 300)
ESSENTIAL_FRAGMENT_USAGE="no"  // "no" for SR data; "yes" for PE data to make bigWig tracks with reconstituted fragments
ESSENTIAL_BAMCOVERAGE="--binSize 1 --skipNonCoveredRegions --normalizeUsing CPM" // NO smoothing should be done for RNAseq

// Alignment with Cellranger (used for main 10X pipeline branch)
ESSENTIAL_TENX_TRANSCRIPTOME="/fsimb/common/genomes/homo_sapiens/10X/grch38_ensembl98/refdata-gex-GRCh38-2020-A"
ESSENTIAL_TENX_REFERENCE="/fsimb/common/genomes/homo_sapiens/10X/grch38_ensembl98/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"
ESSENTIAL_TENX_EXPECTED_CELLS=7000 // cellranger default 
ESSENTIAL_TENX_FASTQDIR=ESSENTIAL_PROJECT + "/rawdata"
ESSENTIAL_TENX_AGGRCSV="aggregation.csv"
ESSENTIAL_TENX_NORMALIZED="mapped"

// Genome annotation
ESSENTIAL_GENESGTF="/fsimb/common/genomes/homo_sapiens/10X/grch38_ensembl98/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
ESSENTIAL_GENESBED="/fsimb/common/genomes/homo_sapiens/10X/grch38_ensembl98/refdata-gex-GRCh38-2020-A/genes/genes.bed"
ESSENTIAL_CHROMSIZES="/fsimb/common/genomes/homo_sapiens/10X/grch38_ensembl98/refdata-gex-GRCh38-2020-A/star/chrNameLength.txt"  // chromosome sizes file
ESSENTIAL_FEATURETYPE="gene_type" //gencode (also 10X Genomics) uses gene_type; ensemble uses gene_biotype
ESSENTIAL_ORG="human"           // UCSC organism
ESSENTIAL_DB="hg38"              // UCSC assembly version
ESSENTIAL_MTGENES=""  // list with gene_id of mitochondrial genes (give path within ESSENTIAL_PROJECT)


// STARR mRNA amplicon library settings
ESSENTIAL_STARR_STAR_REF="/path/to/starrseq_constructs/star_index/"
ESSENTIAL_STARR_GTF="/path/to/starrseq_constructs/constructs.gtf"
ESSENTIAL_STARR_BED="/path/to/starrseq_constructs/constructs.bed"                       // (not currently used)
ESSENTIAL_STARR_CHROMSIZES="/path/to/starrseq_constructs/star_index/chrNameLength.txt"  // chromosome sizes file (not currently used)
//STARR construct polyadenylation signal 3' end sequence (reverse complement) for demultiplexing
// SV40 PAS 3'-UTR final 20 nt used here (actually from -22 to -2, because it ends in 'CA', which would still match the 10X poly-dT primer)
// For reference: SV40 PAS 3'-UTR final 25 nt coding strand sequence: AGCTGCAATAAACAAGTTAACAACA
// For reference: SV40 cleavage site obtained from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC85533/
ESSENTIAL_STARR_FILTERSEQ='TTGTTAACTTGTTTATTGCA'
// Barcodes (used for UMI and barcode extraction in the STARR mRNA processing)
ESSENTIAL_BCPATTERN="CCCCCCCCCCCCCCCCNNNNNNNNNNNN" //10X GEX Read 1 barcode pattern (C = Cell barcode, N = UMI).
ESSENTIAL_WHITELIST="" //the barcode list should be a list of valid barcodes separated by newline (optional)


// Reference genomes for FastqScreen
ESSENTIAL_FASTQSCREEN_PERC=1    // contaminant filter, if a contaminant is consuming at least this percentage of reads in at least one sample, contaminant will be shown in report
ESSENTIAL_FASTQSCREEN_GENOME="Human::/fsimb/common/genomes/homo_sapiens/gencode/release-25_GRCh38.p7/full/index/bowtie2/GRCh38.p7.genome"  //bowtie2 reference for the genome the samples are from, this is used for the fastqscreen
ESSENTIAL_FASTQSCREEN=ESSENTIAL_FASTQSCREEN_GENOME + ",PHIX::/fsimb/common/genomes/phix/19930428/NCBI/index/bowtie2/2.3.4.3/ncbi_phix,ERCC::/fsimb/common/genomes/ERCC/index/bowtie2/2.3.4.3/ERCC92,rRNA::/fsimb/common/genomes/contaminants/fastqscreen_references/rrna/v1/index/bowtie2/2.3.4.3/hs_mm_ce_dm_rn_dr_xt_rRNA,Mycoplasma::/fsimb/common/genomes/contaminants/fastqscreen_references/mycoplasma/v1/index/bowtie2/2.3.4.3/mycoplasma_all_ref,E.coli::/fsimb/common/genomes/Escherichia_coli/ensembl/full/index/bowtie2/Escherichia_coli_str_k_12_substr_dh10b.ASM1942v1.31.dna.genome,B.taurus::/fsimb/common/genomes/bos_taurus/ensembl/3.1/full/index/bowtie2/2.2.9/UMD3.1" //references for fastqscreen to use if run, this are our standard references please include yours 


//Adapter trimming
RUN_CUTADAPT=false
//ESSENTIAL_ADAPTER_SEQUENCE="TruSeqLTHT=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" 
ESSENTIAL_ADAPTER_SEQUENCE="Nextera=CTGTCTCTTATACACATCT" //standard sequence to trim illumina reads
ESSENTIAL_MINADAPTEROVERLAP=10
ESSENTIAL_MINREADLENGTH=30
ESSENTIAL_BASEQUALCUTOFF=20    // trim low-quality ends from reads (if nextseqtrim is true, qualities of terminal G bases are ignored)  
ESSENTIAL_NEXTSEQTRIM=true     // accounts for terminal G bases during base quality trimming incorporated by faulty dark cycles observed with two-color chemistry (as in NextSeq) 

//global vars
PROJECT=ESSENTIAL_PROJECT
LOGS=PROJECT + "/logs"
MAPPED=PROJECT + "/mapped"
TRIMMED=PROJECT + "/trimmed"
SPLITMRNA=PROJECT + "/splitmrna"
QC=PROJECT + "/qc"
REPORTS=PROJECT + "/reports"
RESULTS=PROJECT + "/results"
TMP=PROJECT + "/tmp"
TRACKS=PROJECT + "/tracks"
FUSION=PROJECT + "/fusion"

// optional pipeline stages to include
RUN_TRACKHUB=false
RUN_IN_PAIRED_END_MODE=(ESSENTIAL_PAIRED == "yes")


