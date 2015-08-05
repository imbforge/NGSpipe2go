//vars for task DE_DESeq2 from catalog RNAseq, version 1
DE_DESeq2_TARGETS="targets.txt"     //targets file. Check the bin directory for the format
DE_DESeq2_CONTRASTS="contrasts.txt" //contrasts file. Check the bin directory for the format
DE_DESeq2_MMATRIX="~group"             //simplified formula, MUST contain the intercept
DE_DESeq2_FILTER="NOTUSED"             //whether to filter invariant genes. NOT USED
DE_DESeq2_PREFIX=""       //prefix to be removed from file names if a sample column is not provided in targets.txt
DE_DESeq2_SUFFIX=""       //suffix to be removed from file names if a sample column is not provided in targets.txt
DE_DESeq2_CWD=RESULTS + "/subread-count"             //current working directory
DE_DESeq2_OUTDIR=RESULTS + "/DE_DESeq2" //output filename base pattern

