//vars for task DE_edgeR from catalog RNAseq, version 1
DE_edgeR_TARGETS="targets=targets.txt"       //targets file. Check the bin directory for the format
DE_edgeR_CONTRASTS="contrasts=contrasts.txt" //contrasts file. Check the bin directory for the format
DE_edgeR_MMATRIX="mmatrix=~0+group"          //glm formula with the contrasts to be tested. Check edgeR man
DE_edgeR_FILTER="filter=TRUE"                //whether to filter invariant genes
DE_edgeR_PREFIX="prefix=" //prefix to be removed from file names if a sample column is not provided in targets.txt
DE_edgeR_SUFFIX="suffix=" //suffix to be removed from file names if a sample column is not provided in targets.txt
DE_edgeR_CWD="cwd=" + RESULTS + "/subread-count" //current working directory
DE_edgeR_ROBUST="robust=FALSE"           //robust estimation of the negative binomial dispersion
DE_edgeR_OUTDIR="out=" + RESULTS + "/DE_edgeR"   //output filename base pattern
DE_edgeR_GTF="gtf=" + ESSENTIAL_GENESGTF
DE_edgeR_EXTRA=""   // extra parms to sent to the tool
