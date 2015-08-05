//vars for task RNAtypes from catalog RNAseq, version 1
RNAtypes_FOLDER=MAPPED		// folder containing the mapped files
RNAtypes_PATTERN="duprm\\\\.bam\$"// process all the .bam files
RNAtypes_GTF=ESSENTIAL_GENESGTF2 // the gencode annotation GTF (can be compressed)
RNAtypes_OUT="RNAtypes"		// output file name
RNAtypes_PRE=ESSENTIAL_SAMPLE_PREFIX // prefix to be removed
RNAtypes_SUF=".bam"			// sufix to be removed
RNAtypes_PAIRED=ESSENTIAL_PAIRED		// paired end yes|no
RNAtypes_STRANDED=ESSENTIAL_STRANDED	// strandness yes|no|reverse
RNAtypes_MULTIMAP="NONE"	// ALL|NONE|RANDOM
RNAtypes_FTYPE="exon"		// feature type on GTF file: exon|gene|...
RNAtypes_FTYPECOL="gene_type" // column name  on GTF file containing the biotype information
RNAtypes_CORES=4			// number of cores to use

