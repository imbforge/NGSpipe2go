//vars for task DE_edgeR from catalog RNAseq, version 1
DE_edgeR_TARGETS="targets.txt"     //targets file. Check the bin directory for the format
DE_edgeR_CONTRASTS="contrasts.txt" //contrasts file. Check the bin directory for the format
DE_edgeR_MMATRIX="~0+group"          //glm formula with the contrasts to be tested. Check edgeR man
DE_edgeR_FILTER="TRUE"             //whether to filter invariant genes
DE_edgeR_PREFIX="" //prefix to be removed from file names if a sample column is not provided in targets.txt
DE_edgeR_SUFFIX="" //suffix to be removed from file names if a sample column is not provided in targets.txt
DE_edgeR_CWD=RESULTS + "/subread-count"        //current working directory
DE_edgeR_ROBUST="FALSE"            //robust estimation of the negative binomial dispersion
DE_edgeR_OUTDIR=RESULTS + "/DE_edgeR" //output filename base pattern
DE_edgeR_GTF=ESSENTIAL_GENESGTF

