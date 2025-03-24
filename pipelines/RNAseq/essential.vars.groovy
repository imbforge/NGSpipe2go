// RNA-Seq ESSENTIAL VARIABLES

// Define essential variables here.
// Further module-specific variables can be adjusted in the corresponding ".header" files for each module.

// General parameters
ESSENTIAL_PROJECT="/project/"   // full project directory path
ESSENTIAL_THREADS=4             // number of threads for parallel tasks
ESSENTIAL_SAMPLE_PREFIX=""      // sample name prefix to be trimmed in the results and reports

// Mapping and annotation parameters
ESSENTIAL_STAR_REF="..../star/2.7.3a"     // directory containing all STAR index files
ESSENTIAL_GENESGTF="..../annotation.gtf"  // properly formatted gene annotation GTF file having *different* gene and trancript IDs (e.g. from Ensembl or Gencode); If you use Gencode here, please remember to edit the gtf module variables in the rmats.header file
ESSENTIAL_GENESBED="..../annotation.bed"  // gene annotation BED file 
ESSENTIAL_CHROMSIZES="..../genome.fa.fai" // chromosome sizes file
ESSENTIAL_FEATURETYPE="gene_type"         // Gencode uses "gene_type", Ensembl uses "gene_biotype"
ESSENTIAL_PAIRED="no"           // single-end ("no") or paired-end ("yes") data
ESSENTIAL_READLENGTH=50         // read length - used for optimal mapping over splice-junctions
ESSENTIAL_STRANDED="reverse"    // library prep protocol strandness: no|yes|reverse
ESSENTIAL_ORG="mouse"           // UCSC organism name (human, mouse, fly, worm, yeast, zebrafish)
ESSENTIAL_DB="mm9"              // UCSC assembly version
ESSENTIAL_FILTER_CHR=""         // chromosomes to include in post-mapping analysis
ESSENTIAL_BAMCOVERAGE="--binSize 1 --skipNonCoveredRegions --normalizeUsing CPM"  // size of bins (in bp) and normalisation method for the bigWig tracks

// DESeq2 specific parameters
ESSENTIAL_DESEQ2_FDR=0.01       // FDR significance cutoff in the DESeq2 model (may be increased for noisy or underpowered datasets)
ESSENTIAL_DESEQ2_FC=1           // optional Fold-Change cutoff to incorporate into the DESeq2 model (default "FC=1" means no Fold-Change filter is used)
// Using FC threshold in the DESeq2 model is usually more conservative than post-hoc gene filtering by FC (which should anyway be avoided, see https://doi.org/10.1093/bib/bbab053) 

// Optional pipeline stages to include
RUN_TRACKHUB=false              // prepare a Track Hub for the UCSC genome browser
RUN_FASTQSCREEN=true            // check for contaminations using FastQ Screen
RUN_CUTADAPT=false              // optional read trimming with Cutadapt e.g. if using high read length
RUN_IN_PAIRED_END_MODE=(ESSENTIAL_PAIRED == "yes") // no need to change this line
RUN_RMATS=true // check for differential splicing events using the provided contrasts; mmatrix of contrasts file will be ignored
// In a typical setup, batch effect correction is not necessary for rMATS, as it calculates a ratio of splicing events with in each sample separately and then compares it with the respective replicates. Before comparing with replicates, they also model the variation among the replicates by random effects in a mixed model. In general, the variation accounting technical differences among replicates seems to be taken care by it. 

// FASTQ-Screen parameters
ESSENTIAL_FASTQSCREEN_PERC=1    // contaminant filter, if a contaminant is consuming at least this percentage of reads in at least one sample, contaminant will be shown in report
ESSENTIAL_FASTQSCREEN_GENOME="Organism::Bowtie2_index"  // bowtie2 reference index files (base name) for the genome the samples are from
ESSENTIAL_FASTQSCREEN=ESSENTIAL_FASTQSCREEN_GENOME + ",PHIX::/fsimb/common/genomes/phix/19930428/NCBI/index/bowtie2/2.3.4.3/ncbi_phix,ERCC::/fsimb/common/genomes/ERCC/index/bowtie2/2.3.4.3/ERCC92,rRNA::/fsimb/common/genomes/contaminants/fastqscreen_references/rrna/v1/index/bowtie2/2.3.4.3/hs_mm_ce_dm_rn_dr_xt_rRNA,Mycoplasma::/fsimb/common/genomes/contaminants/fastqscreen_references/mycoplasma/v1/index/bowtie2/2.3.4.3/mycoplasma_all_ref,E.coli::/fsimb/common/genomes/Escherichia_coli/ensembl/full/index/bowtie2/Escherichia_coli_str_k_12_substr_dh10b.ASM1942v1.31.dna.genome,B.taurus::/fsimb/common/genomes/bos_taurus/ensembl/3.1/full/index/bowtie2/2.2.9/UMD3.1"  // bowtie2 reference index files (base name) for common sample contaminants, as well as for PhiX and ERCC library spike-ins

// Adapter trimming with Cutadapt (additional adapter sequences for R1 and/or R2 can be specified in the cutadapt.header file)
// Cutadapt recommends using full length adapter sequences since adapter fragments might occur in the genome
// The sequences to trim for specific libraries can be found here https://knowledge.illumina.com/library-preparation/general/library-preparation-general-reference_material-list/000001314
// The one used most frequently is the Nexterra Transposase Adapter CTGTCTCTTATACACATCT
ESSENTIAL_ADAPTER_SEQUENCE="Illumina=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" // 3'adapter sequence for R1. 
ESSENTIAL_ADAPTER_SEQUENCE_R2="" // in case of paired end sequencing the R2 adapter sequence (which will be trimmed from the 3' end of R2, -A argument in Cutadapt), by default the same sequence as defined in ESSENTIAL_ADAPTER_SEQUENCE is used
ESSENTIAL_MINADAPTEROVERLAP=5
ESSENTIAL_MINREADLENGTH=30
ESSENTIAL_BASEQUALCUTOFF=20  // trim low-quality ends from reads (if nextseqtrim is true, qualities of terminal G bases are ignored)  
ESSENTIAL_NEXTSEQTRIM=true   // accounts for terminal G bases during base quality trimming incorporated by faulty dark cycles observed with two-color chemistry (as in NextSeq) 


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
