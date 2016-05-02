RNAtypes_FOLDER="folder=" + MAPPED	       // folder containing the mapped files
RNAtypes_PATTERN="pattern=duprm\\\\.bam\$" // process all the .bam files
RNAtypes_GTF="gtf=" + ESSENTIAL_GENESGTF2  // the gencode annotation GTF (can be compressed)
<<<<<<< HEAD
RNAtypes_OUT="RNAtypes"		            // output file name prefix the output dir is defined in the module
=======
RNAtypes_OUTNAME="RNAtypes"                // output file name
RNAtypes_OUT="out=" + RNAtypes_OUTNAME     // output file name (parm)
>>>>>>> 09da5224d240850b96c336eca90e9732bac000b9
RNAtypes_PRE="pre=" + ESSENTIAL_SAMPLE_PREFIX // prefix to be removed
RNAtypes_SUF="suf=.bam"			            // sufix to be removed
RNAtypes_PAIRED="paired=" + ESSENTIAL_PAIRED   // paired end yes|no
RNAtypes_STRANDED="stranded=" + ESSENTIAL_STRANDED // strandness yes|no|reverse
RNAtypes_MULTIMAP="multimap=NONE"	        // ALL|NONE|RANDOM
RNAtypes_FTYPE="ftype=exon"		            // feature type on GTF file: exon|gene|...
RNAtypes_FTYPECOL="ftypecol=gene_type"      // column name  on GTF file containing the biotype information
RNAtypes_CORES="cores=" + Integer.toString(ESSENTIAL_THREADS) // number of cores to use
RNAtypes_EXTRA=""

