//vars for task DE_DESeq2 from catalog RNAseq, version 1
DE_DESeq2_TARGETS="targets=targets.txt"       //targets file. Check the bin directory for the format
DE_DESeq2_CONTRASTS="contrasts=contrasts.txt" //contrasts file. Check the bin directory for the format
DE_DESeq2_MMATRIX="mmatrix=~condition"            //simplified formula, MUST contain the intercept
DE_DESeq2_FILTER="filter=TRUE"             //whether to filter invariant genes. NOT USED
DE_DESeq2_PREFIX="prefix="       //prefix to be removed from file names if a sample column is not provided in targets.txt
DE_DESeq2_SUFFIX="suffix="       //suffix to be removed from file names if a sample column is not provided in targets.txt
DE_DESeq2_CWD="cwd=" + RESULTS + "/subread-count"  //current working directory
DE_DESeq2_OUTDIR="out=" + RESULTS + "/DE_DESeq2"   //output filename base pattern. If you change it here, change it also in the module file
DE_DESeq2_EXTRA=""                                 //extra parms to sent to the tool

