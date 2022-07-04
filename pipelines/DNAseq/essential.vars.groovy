// DNA-Seq ESSENTIAL VARIABLES

// Define essential variables here.
// Further module-specific variables can be adjusted in the corresponding ".header" files for each module.

// General parameters
ESSENTIAL_PROJECT="your/project/folder"
ESSENTIAL_SAMPLE_PREFIX="" 
ESSENTIAL_CALL_REGION="" // limit variant calling to a limited region (required), e.g. exonic sequences or chr3 only. Will save a lot of compute time. List of chromosome names is also fine.
ESSENTIAL_THREADS=4
ESSENTIAL_BAMCOVERAGE="--binSize 10 --normalizeUsing CPM"  // deepTools options for making normalised bigWig tracks

// Mapping parameters
ESSENTIAL_BWA_REF="/fsimb/common/genomes/homo_sapiens/gencode/release-35_GRCh38.p13/canonical/index/bwa/0.7.15/GRCh38.canonical_assembly.genome.fa"
ESSENTIAL_PAIRED="yes"        // paired end design
ESSENTIAL_DEDUPLICATION="yes" // remove duplicated reads (default "yes", set "no" only for special designs like amplicon-seq)

// Variant annotation and filtering
ESSENTIAL_FILTERVARIANTS="VQSR" // either "VQSR", "hard-filter" or empty (VQSR is recommended for variant filtering if possible)
ESSENTIAL_KNOWN_VARIANTS="/fsimb/common/genomes/homo_sapiens/gatk_bundle_hg38_v0/dbsnp_146.hg38.vcf.gz" // KnownVariants could be dbSNP from GATK resource bundle. This is crucial for BaseQualityRecalibration step.
ESSENTIAL_HAPMAP_VARIANTS="/fsimb/common/genomes/homo_sapiens/gatk_bundle_hg38_v0/hapmap_3.3.hg38.vcf.gz" // varirants provided by the GATK bundle. Essential for Variant Score Recalibration
ESSENTIAL_OMNI_VARIANTS="/fsimb/common/genomes/homo_sapiens/gatk_bundle_hg38_v0/1000G_omni2.5.hg38.vcf.gz" // varirants provided by the GATK bundle. Essential for Variant Score Recalibration
ESSENTIAL_MILLS_VARIANTS="/fsimb/common/genomes/homo_sapiens/gatk_bundle_hg38_v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" // varirants provided by the GATK bundle. Essential for Variant Score Recalibration
ESSENTIAL_THOUSAND_GENOMES_VARIANTS="/fsimb/common/genomes/homo_sapiens/gatk_bundle_hg38_v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz" // varirants provided by the GATK bundle. Essential for Variant Score Recalibration
ESSENTIAL_SNPEFF_GENOME="GRCh38.p13.RefSeq" // snpEff resource to be used for functional variant annotation

// FASTQ-Screen parameters (optional)
ESSENTIAL_FASTQSCREEN_PERC=1    // contaminant filter, if a contaminant is consuming at least this percentage of reads in at least one sample, contaminant will be shown in report
ESSENTIAL_FASTQSCREEN_GENOME="Human::/fsimb/common/genomes/homo_sapiens/gencode/release-35_GRCh38.p13/full/index/bowtie2/2.3.4.3/GRCh38.p13"  //bowtie2 reference for the genome the samples are from, this is used for the fastqscreen
ESSENTIAL_FASTQSCREEN=ESSENTIAL_FASTQSCREEN_GENOME + ",PHIX::/fsimb/common/genomes/phix/19930428/NCBI/index/bowtie2/2.3.4.3/ncbi_phix,ERCC::/fsimb/common/genomes/ERCC/index/bowtie2/2.3.4.3/ERCC92,rRNA::/fsimb/common/genomes/contaminants/fastqscreen_references/rrna/v1/index/bowtie2/2.3.4.3/hs_mm_ce_dm_rn_dr_xt_rRNA,Mycoplasma::/fsimb/common/genomes/contaminants/fastqscreen_references/mycoplasma/v1/index/bowtie2/2.3.4.3/mycoplasma_all_ref,E.coli::/fsimb/common/genomes/Escherichia_coli/ensembl/full/index/bowtie2/Escherichia_coli_str_k_12_substr_dh10b.ASM1942v1.31.dna.genome,B.taurus::/fsimb/common/genomes/bos_taurus/ensembl/3.1/full/index/bowtie2/2.2.9/UMD3.1" //references for fastqscreen to use if run, this are our standard references please include yours 

// Adapter trimming with Cutadapt (optional; additional adapter sequences for R1 and/or R2 can be specified in the cutadapt.header file)
ESSENTIAL_ADAPTER_SEQUENCE="AmpliSeq=CTGTCTCTTATACACATCT" // 3'adapter sequence for R1. 
ESSENTIAL_MINADAPTEROVERLAP=5
ESSENTIAL_MINREADLENGTH=30
ESSENTIAL_BASEQUALCUTOFF=20  // trim low-quality ends from reads (if nextseqtrim is true, qualities of terminal G bases are ignored)  
ESSENTIAL_NEXTSEQTRIM=true   // accounts for terminal G bases during base quality trimming incorporated by faulty dark cycles observed with two-color chemistry (as in NextSeq) 


// Optional pipeline stages to include
RUN_TRACKHUB=false              // prepare a Track Hub for the UCSC genome browser
RUN_FASTQSCREEN=true            // check for contaminations using FastQ Screen
RUN_CUTADAPT=false              // optional read trimming with Cutadapt e.g. if using high read length
RUN_SNPEFF=true
RUN_IN_PAIRED_END_MODE=(ESSENTIAL_PAIRED == "yes") // no need to change this line
RUN_RMDUPS=(ESSENTIAL_DEDUPLICATION == "yes")

//project folders
PROJECT=ESSENTIAL_PROJECT
LOGS=PROJECT + "/logs"
TRIMMED=PROJECT + "/trimmed"
MAPPED=PROJECT + "/mapped"
QC=PROJECT + "/qc"
REPORTS=PROJECT + "/reports"
RESULTS=PROJECT + "/results"
TMP=PROJECT + "/tmp"
TRACKS=PROJECT + "/tracks"


