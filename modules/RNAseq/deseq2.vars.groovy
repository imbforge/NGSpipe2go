DE_DESeq2_TARGETS="targets=targets.txt"       //targets file. Check the bin directory for the format
DE_DESeq2_CONTRASTS="contrasts=contrasts.txt" //contrasts file. Check the bin directory for the format
DE_DESeq2_MMATRIX="mmatrix=~group"        //formula for the linear model
DE_DESeq2_FILTER="filter=TRUE"                //whether to filter invariant genes
DE_DESeq2_PREFIX="prefix="                    //prefix to be removed from file names
DE_DESeq2_SUFFIX="suffix="                    //suffix to be removed from file names
DE_DESeq2_CWD="cwd=" + RESULTS + "/subread-count"  //current working directory
DE_DESeq2_OUTDIR="out=" + RESULTS + "/DE_DESeq2"   //output filename base pattern. If you change it here, change it also in the module file
DE_DESeq2_GENES=ESSENTIAL_GENESGTF		   //adding a parameter for the gtf file (DESeq2 module)
DE_DESeq2_BASE="base="				   //base level log2foldchanges should be calulated for e.g. wt or control has to be a group id within the targets file
DE_DESeq2_EXTRA=""                                 //extra parms to sent to the tool

