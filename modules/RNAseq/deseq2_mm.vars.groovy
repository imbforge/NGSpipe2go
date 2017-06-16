DE_DESeq2_MM_TARGETS=DE_DESeq2_TARGETS       //targets file. Check the bin directory for the format
DE_DESeq2_MM_CONTRASTS=DE_DESeq2_CONTRASTS   //contrasts file. Check the bin directory for the format
DE_DESeq2_MM_MMATRIX=DE_DESeq2_MMATRIX       //formula for the linear model
DE_DESeq2_MM_FILTER=DE_DESeq2_FILTER         //whether to perform automatic independent filtering of lowly expressed genes
DE_DESeq2_MM_PREFIX=DE_DESeq2_PREFIX         //prefix to be removed from file names
DE_DESeq2_MM_SUFFIX=DE_DESeq2_SUFFIX         //suffix to be removed from file names
DE_DESeq2_MM_DUPRADAR_OUTDIR=DUPRADAR_OUTDIR.replaceFirst("outdir=", "") //where the dupradar left the output tables
DE_DESeq2_MM_CWD="cwd=" + RESULTS + "/subread-count_MM"  //directory where the .tsv files with the gene counts are located
DE_DESeq2_MM_OUTDIR="out=" + RESULTS + "/DE_DESeq2_MM"   //output filename base pattern. If you change it here, change it also in the module file
DE_DESeq2_MM_GENES=DE_DESeq2_GENES           //gtf file to calculate RPKM
DE_DESeq2_MM_BASE=DE_DESeq2_BASE             //base level log2foldchanges should be calulated for e.g. wt or control has to be a group id within the targets file
DE_DESeq2_MM_EXTRA=DE_DESeq2_EXTRA           //extra parms to sent to the tool

