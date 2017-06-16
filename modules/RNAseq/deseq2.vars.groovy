DE_DESeq2_TARGETS="targets=targets.txt"       //targets file. Check the bin directory for the format
DE_DESeq2_CONTRASTS="contrasts=contrasts.txt" //contrasts file. Check the bin directory for the format
DE_DESeq2_MMATRIX="mmatrix=~group"        //formula for the linear model
DE_DESeq2_FILTER="filter=TRUE"                //whether to perform automatic independent filtering of lowly expressed genes
DE_DESeq2_PREFIX="prefix="                    //prefix to be removed from file names
DE_DESeq2_SUFFIX="suffix="                    //suffix to be removed from file names
DE_DESeq2_CWD="cwd=" + RESULTS + "/subread-count"  //directory where the .tsv files with the gene counts are located
DE_DESeq2_OUTDIR="out=" + RESULTS + "/DE_DESeq2"   //output filename base pattern. If you change it here, change it also in the module file
DE_DESeq2_GENES="gtf=" + ESSENTIAL_GENESGTF		      //gtf file to calculate RPKM
DE_DESeq2_BASE="base="				   //base level log2foldchanges should be calulated for e.g. wt or control has to be a group id within the targets file
DE_DESeq2_EXTRA=""                                 //extra parms to sent to the tool

