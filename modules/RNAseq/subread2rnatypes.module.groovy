//rule for task feature_count from Subread package, version 1
//desc: Counting reads in features with featureCounts

// example content of env.sh
// export PATH=${TOOL_DEPENDENCIES}/subread/1.5.0/bin:$PATH

rnatypes = {
	doc title: "subread_count",
	desc:  "Counting gene biotypes in features with feature-count out of the subread package",
	constraints: """Default: strand specific counting.""",
	bpipe_version: "tested with bpipe 0.9.8.7",
	author: "Oliver Drechsel"
	
	output.dir  = RNATYPES_OUTDIR
	def RNATYPES_FLAGS = "-F GTF " +
						RNATYPES_GENESGTF + " " +
						RNATYPES_FEATURE + " " +
						RNATYPES_ACCUMULATE + " " +
					    RNATYPES_CORES    + " " +
                        RNATYPES_EXTRA    + " "
						
	
	if(RNATYPES_PAIRED == "yes") {
		RNATYPES_FLAGS = "-p " + RNATYPES_FLAGS
	}
	
	// no|yes|reverse
	if(RNATYPES_STRANDED == "no") {
		RNATYPES_FLAGS = "-s 0 " + RNATYPES_FLAGS
	}
	else if (RNATYPES_STRANDED == "yes") {
		RNATYPES_FLAGS = "-s 1 " + RNATYPES_FLAGS
	}
	else {
		RNATYPES_FLAGS = "-s 2 " + RNATYPES_FLAGS
	}
	
	// run the chunk
	transform(".bam") to ("_readcounts.tsv") {
		exec """
			module load subread &&
			
			echo 'VERSION INFO'  1>&2 ;
			echo \$(featureCounts 2>&1 | grep Version | cut -d' ' -f2) 1>&2 ;
			echo '/VERSION INFO' 1>&2 ;
			
			featureCounts $RNATYPES_FLAGS -o ${output}_tmp $input 2> ${output.prefix}_rnatypeslog.stderr;
			cut -f1,6,7 ${output}_tmp > $output;
			rm ${output}_tmp;	

		""","rnatypes"
	}
}
