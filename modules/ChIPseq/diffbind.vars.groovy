DIFFBIND_TARGETS="targets=targets_diffbind.txt"        //targets file. Check the bin directory for the format
DIFFBIND_CONTRASTS="contrasts=contrasts_diffbind.txt" //contrasts file. Check the bin directory for the format
DIFFBIND_CWD="cwd=" + PROJECT                //current working directory
DIFFBIND_OUTDIR="out=" + RESULTS + "/diffbind" //output directory
DIFFBIND_BAMS="bams=" + MAPPED               // directory with the bam files
DIFFBIND_FRAGSIZE="fragsize=" + Integer.toString(ESSENTIAL_FRAGLEN) // average fragment size
DIFFBIND_ANNOTATE="annotate=TRUE"            // annotate peaks after differential bindign analysis?
DIFFBIND_TSS="tss='c(-3000,3000)'"             // region around the TSS to be considered as promoter
DIFFBIND_TXDB="txdb=" + ESSENTIAL_TXDB       // Bioconductor transcript database, for annotation
DIFFBIND_ANNODB="annodb=" + ESSENTIAL_ANNODB // Bioconductor gene annotation database
DIFFBIND_PAIRED="pe=" + ((ESSENTIAL_PAIRED == "yes") ? "TRUE" : "FALSE") // Paired end experiment?
DIFFBIND_EXTRA=""                            //extra parms to sent to the tool
