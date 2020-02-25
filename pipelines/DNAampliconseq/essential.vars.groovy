ESSENTIAL_PROJECT="your_project_folder"
ESSENTIAL_THREADS=4

ESSENTIAL_EXPDESIGN="amplicon1"  // design of DNA Amplicon-Seq experiment.
// the following experimental designs are implemented:
// "amplicon1": Paired end sequencing with overlapping reads to be merged.
//              Amplicon sequence contains cell barcodes to count and UMIs for deduplication.
//              Other elements are optional, but must fit the regular expression given in "ESSENTIAL_BCPATTERN".
// "amplicon2": as "amplicon1", but without UMIs. Remark: the regular expression still needs an UMI segment somewhere
//              for reasons of compatibility with umi_tools, though of length zero: (?P<umi_1>.{0}).
// "amplicon3": Paired end sequencing with non-overlapping reads. No merging of read pairs before barcode extraction.
//              Contains two independent cell barcodes but no UMIs. The occurrences of all CB combinations shall be counted.


// use predefined BC patterns for known experiment design. If other sequences are used define patterns below.
switch(ESSENTIAL_EXPDESIGN) {
    case "amplicon1":
       PREDEF_BCPATTERN="\"(?P<discard_1>.{0,7})(?P<umi_1>.{8})(CTACAACATCCAGAAGGAGTCGACCCTGCACCTGGTCCTGCGTCTGAGAGGTGGT){s<=2}(?P<cell_1>.{6})(GGATCCGGAGCTTGGCTGTTGCCCGTCTCACTGGTGAAAAGAAAAACCACCCTGGC){s<=2}(?P<discard_2>.*)\"";
       PREDEF_BCPATTERN_2="";
       break;
    case "amplicon2":
       PREDEF_BCPATTERN="\"(?P<discard_1>.{0,10})(?P<umi_1>.{0})(AGGAGTCCACCTTACATCTTGTGCTAAGGCTAAGAGGTGGT){s<=2}(?P<cell_1>.{6})(GGATCCGGAGCTTGGCTGTTGCCCGTCTCACTGGTGAAAAGAAAAACCACCCTGGCGCCCAATA){s<=2}(?P<discard_2>.*)\"";
       PREDEF_BCPATTERN_2="";
       break;
    case "amplicon3":
       PREDEF_BCPATTERN="\"(?P<discard_1>.*)(?P<umi_1>.{0})(ATGGACGAATTGTACAAGGAACAGAAGTTGATTTCTGAAGAAGACCTCGGTTCTGGATCAGGTTCA){s<=2}(?P<cell_1>.{54,120})(TAACTCCAGGAGTATATAAA){s<=2}(?P<discard_2>.*)\"";
       PREDEF_BCPATTERN_2="\"(?P<discard_1>.*)(?P<umi_1>.{0})(CGAATTCAAGCTTAGATCTGATATCGGTACC){s<=2}(?P<cell_1>.{26})(ATAACTTCGTATAGCATACATTATACGAAGTTAT){s<=2}(?P<discard_2>.*)\"";
       break;
    default:
       PREDEF_BCPATTERN="CCCCCCCNNNNNNNN"; // define number of barcode (C) and UMI (N) bases and set ESSENTIAL_EXTRACTMETHOD="string"
       PREDEF_BCPATTERN_2="";
       break;
}

ESSENTIAL_BCPATTERN=PREDEF_BCPATTERN
ESSENTIAL_BCPATTERN_2=PREDEF_BCPATTERN_2

ESSENTIAL_EXTRACTMETHOD="regex" // either "string" or "regex"
ESSENTIAL_WHITELIST="" //the barcode list should be a list of valid barcodes separated by newline


//global vars that will be reused in some global vars
PROJECT=ESSENTIAL_PROJECT
LOGS=PROJECT + "/logs"
MAPPED=PROJECT + "/mapped"
QC=PROJECT + "/qc"
REPORTS=PROJECT + "/reports"
RESULTS=PROJECT + "/results"
TMP=PROJECT + "/tmp"
TRACKS=PROJECT + "/tracks"

// optional pipeline stages to include
RUN_IN_PAIRED_END_MODE=(ESSENTIAL_EXPDESIGN in ["amplicon1","amplicon2","amplicon3"])
RUN_PEAR=(ESSENTIAL_EXPDESIGN in ["amplicon1","amplicon2"])
RUN_MPSprofiling=(ESSENTIAL_EXPDESIGN in ["amplicon1","amplicon2"])