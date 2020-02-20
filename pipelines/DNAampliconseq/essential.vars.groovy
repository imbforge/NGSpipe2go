ESSENTIAL_PROJECT="/fsimb/groups/imb-bioinfocf/projects/khmelinskii/imb_khmelinskii_2019_02_nieto_amplicon_test"
ESSENTIAL_THREADS=4
ESSENTIAL_PAIRED="yes"

ESSENTIAL_EXPDESIGN="amplicon1"  // design of DNA Amplicon-Seq experiment.
// the following experimental designs are implemented:
// "amplicon1": Paired end sequencing with overlapping reads to be merged.
//              Amplicon sequence contains cell barcodes to count and UMIs for demultiplexing.
//              Other elements are optional, but must fit the regular expression given in "ESSENTIAL_BCPATTERN".
// "amplicon2": as "amplicon1", but with no UMIs. Remark: the regular expression still needs an UMI segment anywhere
//              for reasons of compatibility with umi_tools, though of length zero: (?P<umi_1>.{0}).
// "amplicon3": Paired end sequencing with non-overlapping reads. No merging of read pairs before barcode extraction.
//              Contains two independent cell barcodes but no UMIs. The occurrences of all CB combinations shall be counted.
ESSENTIAL_EXTRACTMETHOD="regex" // either "string" or "regex"
ESSENTIAL_WHITELIST="" //the barcode list should be a list of valid barcodes separated by newline
// ESSENTIAL_BCPATTERN="CCCCCCCNNNNNNNN" //barcode pattern as it is present in MARS-Seq data.
ESSENTIAL_BCPATTERN="\"(?P<discard_1>.{0,7})(?P<umi_1>.{8})(CTACAACATCCAGAAGGAGTCGACCCTGCACCTGGTCCTGCGTCTGAGAGGTGGT){s<=2}(?P<cell_1>.{6})(GGATCCGGAGCTTGGCTGTTGCCCGTCTCACTGGTGAAAAGAAAAACCACCCTGGC){s<=2}(?P<discard_2>.*)\"" //amplicon1
//ESSENTIAL_BCPATTERN="\"(?P<discard_1>.{0,20})(?P<umi_1>.{0})(CTACAACATCCAGAAGGAGTCGACCCTGCACCTGGTCCTGCGTCTGAGAGGTGGT){s<=2}(?P<cell_1>.{6})(GGATCCGGAGCTTGGCTGTTGCCCGTCTCACTGGTGAAAAGAAAAACCACCCTGGC){s<=2}(?P<discard_2>.*)\"" //amplicon2
// ESSENTIAL_BCPATTERN="\"(?P<discard_1>.{0,10})(?P<umi_1>.{0})(AGGAGTCCACCTTACATCTTGTGCTAAGGCTAAGAGGTGGT){s<=2}(?P<cell_1>.{6})(GGATCCGGAGCTTGGCTGTTGCCCGTCTCACTGGTGAAAAGAAAAACCACCCTGGCGCCCAATA){s<=2}(?P<discard_2>.*)\""  // amplicon3


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
