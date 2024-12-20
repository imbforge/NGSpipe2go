// smallRNA-Seq ESSENTIAL VARIABLES 

// Define essential variables here. 
// Further module-specific variables can be adjusted in the corresponding ".header" files for each module. 

// General parameters
ESSENTIAL_PROJECT="/project/"
ESSENTIAL_THREADS=4
ESSENTIAL_SAMPLE_PREFIX=""
ESSENTIAL_SMALLRNA="miRNA"     // type of smallRNA to be analyzed (types of smallRNA with only few instances should not be analyzed separately
                               //                                  since no proper normalization might be possible)

// Mapping and annotation parameters
ESSENTIAL_BOWTIE_REF="..../bowtie/1.1.2/..."
ESSENTIAL_GENOME_REF="..../genome.fa" // necessary for miRDeep2, not currently used

ESSENTIAL_GENESGTF="..../annotation.gtf"
ESSENTIAL_MIRNAGFF="..../annotation/miRBase/....gff3"
ESSENTIAL_FEATURETYPE="gene_type" //gencode uses gene_type; ensemble uses gene_biotype
ESSENTIAL_PAIRED="no"          // paired end design
ESSENTIAL_STRANDED="yes"       // strandness: no|yes|reverse
ESSENTIAL_ORG="mouse"          // organism
ESSENTIAL_BAMCOVERAGE="--binSize 1 --skipNonCoveredRegions --normalizeUsing CPM" // NO smoothing should be done for RNAseq
// size of bins (in bases)  and method to normalize the number of reads per bin to generate bigwig file.

// Parameters describing the reads
ESSENTIAL_READLENGTH="75"             // actual read length in original raw data (incl. insert, UMIs, adapter)
ESSENTIAL_UMI_LENGTH_LEFT="4"         // length (bp) of left UMI (in case of no UMIs, please set TRIM_UMIS=false and do NOT set ESSENTIAL_UMI_LENGTH_LEFT="0")
ESSENTIAL_UMI_LENGTH_RIGHT="4"        // length (bp) of right UMI (in case of no UMIs, please set TRIM_UMIS=false and do NOT set ESSENTIAL_UMI_LENGTH_RIGHT="0")
ESSENTIAL_MINREADLENGTH_EXCL_UMI="18" // remaining read length excl UMIs
ESSENTIAL_MINREADLENGTH=ESSENTIAL_MINREADLENGTH_EXCL_UMI.toInteger() + ESSENTIAL_UMI_LENGTH_LEFT.toInteger() + ESSENTIAL_UMI_LENGTH_RIGHT.toInteger() // remaining read length plus UMIs
ESSENTIAL_MINADAPTEROVERLAP=5       // minimal overlap with adapter
ESSENTIAL_MAXREADLENGTH=ESSENTIAL_READLENGTH.toInteger() - ESSENTIAL_MINADAPTEROVERLAP.toInteger() // maximal read length to keep (all w/o adapter are discarded)
                                                                                                                     // this parameter needs to changed in case of long smallRNAs, 
                                                                                                                     // which might not reach the adapter or second UMIs
// Cutadapt recommends using full length adapter sequences since adapter fragments might occur in the genome
ESSENTIAL_ADAPTER_SEQUENCE="TruSeqSmallRNA=TGGAATTCTCGGGTGCCAAGG" // needed for cutadapt adapter trimming
ESSENTIAL_ADAPTER_SEQUENCE_R2="" // in case of paired end sequencing the R2 adapter sequence (which will be trimmed from the 3' end of R2, -A argument in Cutadapt), by default the same sequence as defined in ESSENTIAL_ADAPTER_SEQUENCE is used
ESSENTIAL_BASEQUALCUTOFF=0     // trim low-quality ends from reads (if nextseqtrim is true, qualities of terminal G bases are ignored)  
                               // keep this parameter at 0 to switch off low-quality base pair trimming at the ends of reads
                               // (smallRNAseq reads often have UMIs at the end, which should not be removed)
ESSENTIAL_NEXTSEQTRIM=false    // accounts for terminal G bases during base quality trimming incorporated by faulty dark cycles observed with two-color chemistry (as in NextSeq) 
ESSENTIAL_MINIMAL_QUAL="20"    // all reads with any base quality below this value will be removed during quality filtering
                               // (this parameter needs to adjusted in case of long smallRNAs, which can tolerate few bases with low sequence quality)

// DESeq2 specific parameters 
ESSENTIAL_DESEQ2_FDR=0.01    // FDR significance cutoff in the DESeq2 model (may be increased for noisy or underpowered datasets) 
ESSENTIAL_DESEQ2_FC=1        // optional Fold-Change cutoff to incorporate into the DESeq2 model (default "FC=1" means no Fold-Change filter is used) 
// Using FC threshold in the DESeq2 model is usually more conservative than post-hoc gene filtering by FC (which should anyway be avoided, see https://doi.org/10.1093/bib/bbab053) 

// FASTQ-Screen parameters
ESSENTIAL_FASTQSCREEN_PERC="1" // contaminant filter, if a contaminant is consuming at least this percentage of reads in at least one sample, contaminant will be shown in report 
ESSENTIAL_FASTQSCREEN_GENOME="Organism::Bowtie2_index"  // bowtie2 reference index files (base name) for the genome the samples are from
ESSENTIAL_FASTQSCREEN=ESSENTIAL_FASTQSCREEN_GENOME + ",PHIX::/fsimb/common/genomes/phix/19930428/NCBI/index/bowtie2/2.3.4.3/ncbi_phix,ERCC::/fsimb/common/genomes/ERCC/index/bowtie2/2.3.4.3/ERCC92,rRNA::/fsimb/common/genomes/contaminants/fastqscreen_references/rrna/v1/index/bowtie2/2.3.4.3/hs_mm_ce_dm_rn_dr_xt_rRNA,Mycoplasma::/fsimb/common/genomes/contaminants/fastqscreen_references/mycoplasma/v1/index/bowtie2/2.3.4.3/mycoplasma_all_ref,E.coli::/fsimb/common/genomes/Escherichia_coli/ensembl/full/index/bowtie2/Escherichia_coli_str_k_12_substr_dh10b.ASM1942v1.31.dna.genome,B.taurus::/fsimb/common/genomes/bos_taurus/ensembl/3.1/full/index/bowtie2/2.2.9/UMD3.1"  // bowtie2 reference index files (base name) for common sample contaminants, as well as for PhiX and ERCC library spike-ins

// vars for mirDeep2 (module NOT well tested!)
//ESSENTIAL_SPECIES="Mouse"      // necessary for miRDeep2, used to refer to UCSC
//ESSENTIAL_MATURE_MIRNA="..../annotation/miRBase_release21/mmu_mature_oneID.fa"
//ESSENTIAL_HAIRPIN_MIRNA="..../annotation/miRBase_release21/mmu_hairpin_oneID.fa"

// project folders
PROJECT=ESSENTIAL_PROJECT
LOGS=PROJECT + "/logs"
LOGS_MY=PROJECT + "/logs_my"
QC=PROJECT + "/qc"
REPORTS=PROJECT + "/reports"
RESULTS=PROJECT + "/results"
TRIMMED=PROJECT + "/rawdata_processed"
MAPPED=PROJECT + "/mapped"
TMP=PROJECT + "/tmp"
TRACKS=PROJECT + "/tracks"

// optional pipeline stages to include
REMOVE_DUPLICATES=false    // remove duplicate reads
TRIM_UMIS=true             // trim UMIs from read starts and ends (trimming needed when e.g. libraries were prepared using the NEXTflex Small RNA-Seq Kit v3)
RUN_FASTQSCREEN=true       // check for contaminations using FastQ Screen
RUN_MATUREMIRNA_ANALYSIS=true // set to false if no mature miRNAs should be analyzed
RUN_IN_PAIRED_END_MODE=(ESSENTIAL_PAIRED == "yes") // no need to change this line

