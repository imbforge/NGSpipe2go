ESSENTIAL_PROJECT="/your/project/folder/"
ESSENTIAL_THREADS=4
ESSENTIAL_SAMPLE_PREFIX=""

ESSENTIAL_EXPDESIGN="amplicon3"  // design of DNA Amplicon-Seq experiment.
// the following experimental designs are implemented:
// "amplicon1": Paired end sequencing with overlapping reads to be merged.
//              Amplicon sequence contains cell barcodes to count and UMIs for deduplication.
//              Other elements are optional, but must fit the regular expression given in "ESSENTIAL_BCPATTERN".
// "amplicon2": as "amplicon1", but without UMIs. Remark: the regular expression still needs an UMI segment
//              for reasons of compatibility with umi_tools, though of length zero: (?P<umi_1>.{0}).
// "amplicon3": As "amplicon2", but with an additional barcode in read2 for sample demultiplexing (i.e. 2 independent cell barcodes but no UMIs).
//              The 2nd barcode will be extracted by an additional umi_tools extract step to keep it separated from 1st barcode and is copied into 2nd column of count file.
//              (if both barcodes were extracted in a single umi_tools extract step they are merged in the read name and not separated by "_").
//              The occurrences of all CB combinations shall be counted.
// "amplicon4": Paired end sequencing with non-overlapping reads. No pear assembly of read pairs before barcode extraction.
//              Contains two independent cell barcodes but no UMIs. The occurrences of all CB combinations shall be counted.


// use predefined BC patterns for known experiment design. If other sequences are used define patterns below.
switch(ESSENTIAL_EXPDESIGN) {
    case "amplicon1":
       PREDEF_BCPATTERN="\"(?P<discard_1>.{0,7})(?P<umi_1>.{8})(CTACAACATCCAGAAGGAGTCGACCCTGCACCTGGTCCTGCGTCTGAGAGGTGGT){s<=2}(?P<cell_1>.{6})(GGATCCGGAGCTTGGCTGTTGCCCGTCTCACTGGTGAAAAGAAAAACCACCCTGGC){s<=2}(?P<discard_2>.*)\"";
       PREDEF_BCPATTERN_2="";
       break;
    case "amplicon2":
       PREDEF_BCPATTERN="\"(?P<discard_1>.{0,6})(?P<umi_1>.{0})(AGGAGTCCACCTTACATCTTGTGCTAAGGCTAAGAGGTGGT){s<=2}(?P<cell_1>.{6})(GGATCCGGAGCTTGGCTGTTGCCCGTCTCACTGGTGAAAAGAAAAACCACCCTGGCGCCCAATA){s<=2}(?P<discard_2>.*)\"";
       PREDEF_BCPATTERN_2="";
       break;
    case "amplicon3":
       PREDEF_BCPATTERN="\"(?P<discard_1>.{9,10})(?P<umi_1>.{0})(AGGAGTCCACCTTACATCTTGTGCTAAGGCTAAGAGGTGGT){s<=2}(?P<cell_1>.{6})(.*)\"";; // mind to use (.*) in the end to pass through rest of sequence to 2nd call of umi_tools extract
       PREDEF_BCPATTERN_2="\"(?P<discard_1>.{0,50})(?P<umi_1>.{0})(GGATCCGGAGCTTGGCTGTTGCCCGTCTCACTGGTGAAAAGAAAAACCACCCTGGCGCCCAATA){s<=3}(?P<cell_1>.{9,13})\"";
       break;
    case "amplicon4":
       PREDEF_BCPATTERN="\"(?P<discard_1>.{9,10})(?P<umi_1>.{0})(AGGAGTCCACCTTACATCTTGTGCTAAGGCTAAGAGGTGGT){s<=2}(?P<cell_1>.{6})(?P<discard_2>.*)\"";
       PREDEF_BCPATTERN_2="\"(?P<cell_1>.{9,13})(?P<umi_1>.{0})(TATTGGGCGCCAGGGTGGTTTTTCTTTTCACCAGTGAGACGGGCAACAGCCAAGCTCCGGATCC){s<=3}(?P<discard_1>.*)\"";
       break;
    default:
       PREDEF_BCPATTERN="CCCCCCCNNNNNNNN"; // define number of barcode (C) and UMI (N) bases and set ESSENTIAL_EXTRACTMETHOD="string"
       PREDEF_BCPATTERN_2="";
       break;
}

ESSENTIAL_BCPATTERN=PREDEF_BCPATTERN
ESSENTIAL_BCPATTERN_2=PREDEF_BCPATTERN_2

ESSENTIAL_EXTRACTMETHOD="regex" // either "string" or "regex"

// whitelist options
ESSENTIAL_WHITELIST="" //the barcode list should be a list of valid barcodes separated by newline. Barcodes will be filtered (and corrected) against the whitelist.
ESSENTIAL_WHITELIST2="whitelist2.txt" // list of valid barcodes for 2nd run of umi_tools extract (if needed)
ESSENTIAL_EXTRACT_WHITELIST=false  // extract whitelist from data instead providing an external whitelist.
ESSENTIAL_CORRECT_CB=false // correct cell barcodes to whitelist (alternatives must be given in whitelist). If false, just filter against whitelist. Omitted if no whitelist is given (neither external or extracted).
ESSENTIAL_EXTRACT_WHITELIST2=false // // extract whitelist in optional 2nd run of umi_tools extract
ESSENTIAL_CORRECT_CB2=false // correct cell barcodes to whitelist2


// FASTQ-Screen parameters
ESSENTIAL_FASTQSCREEN_PERC=1    // contaminant filter, if a contaminant is consuming at least this percentage of reads in at least one sample, contaminant will be shown in report
ESSENTIAL_FASTQSCREEN_GENOME="Human::/fsimb/common/genomes/homo_sapiens/gencode/release-35_GRCh38.p13/full/index/bowtie2/2.3.4.3/GRCh38.p13"  //bowtie2 reference for the genome the samples are from, this is used for the fastqscreen
ESSENTIAL_FASTQSCREEN=ESSENTIAL_FASTQSCREEN_GENOME + ",Yeast::/fsimb/common/genomes/saccharomyces_cerevisiae/ensembl/R64/canonical/index/bowtie2/2.3.4.3/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel,PHIX::/fsimb/common/genomes/phix/19930428/NCBI/index/bowtie2/2.3.4.3/ncbi_phix,ERCC::/fsimb/common/genomes/ERCC/index/bowtie2/2.3.4.3/ERCC92,rRNA::/fsimb/common/genomes/contaminants/fastqscreen_references/rrna/v1/index/bowtie2/2.3.4.3/hs_mm_ce_dm_rn_dr_xt_rRNA,Mycoplasma::/fsimb/common/genomes/contaminants/fastqscreen_references/mycoplasma/v1/index/bowtie2/2.3.4.3/mycoplasma_all_ref,E.coli::/fsimb/common/genomes/Escherichia_coli/ensembl/full/index/bowtie2/Escherichia_coli_str_k_12_substr_dh10b.ASM1942v1.31.dna.genome,B.taurus::/fsimb/common/genomes/bos_taurus/ensembl/3.1/full/index/bowtie2/2.2.9/UMD3.1" //references for fastqscreen to use if run, this are our standard references please include yours 


//Adapter trimming with Cutadapt (optional).
RUN_CUTADAPT=true
ESSENTIAL_ADAPTER_SEQUENCE="IlluminaR1=AGATCGGAAGAGCACACGTCTGAAC" // standard sequence to trim illumina reads (IlluminaR2=AGATCGGAAGAGCGTCGTGTAGGGA, set in cutadapt header if needed)
ESSENTIAL_MINADAPTEROVERLAP=3  // minimal overlap of the read and the adapter for an adapter to be found (cutadapt default 3)
ESSENTIAL_MINREADLENGTH=15     // minimal length of reads to be kept (cutadapt default 0)
ESSENTIAL_BASEQUALCUTOFF=20  // trim low-quality ends from reads (if nextseqtrim is true, qualities of terminal G bases are ignored)  
ESSENTIAL_NEXTSEQTRIM=true   // accounts for terminal G bases during base quality trimming incorporated by faulty dark cycles observed with two-color chemistry (as in NextSeq) 



//global vars that will be reused in some global vars
PROJECT=ESSENTIAL_PROJECT
LOGS=PROJECT + "/logs"
TRIMMED=PROJECT + "/trimmed"
MAPPED=PROJECT + "/mapped"
QC=PROJECT + "/qc"
REPORTS=PROJECT + "/reports"
RESULTS=PROJECT + "/results"
TMP=PROJECT + "/tmp"
TRACKS=PROJECT + "/tracks"

// optional pipeline stages to include
RUN_IN_PAIRED_END_MODE=(ESSENTIAL_EXPDESIGN in ["amplicon1","amplicon2","amplicon3", "amplicon4"])
RUN_PEAR=(ESSENTIAL_EXPDESIGN in ["amplicon1","amplicon2", "amplicon3"])
RUN_MPSprofiling=(ESSENTIAL_EXPDESIGN in ["amplicon1","amplicon2", "amplicon3"])
