//rule for task RNAtypes from catalog RNAseq, version 1
//desc: Analysis of duplication rate on RNAseq analysis
RNAtypes = {
	doc title: "RNAtypes",
		desc:  "analysis of duplication rate on RNAseq analysis",
		constraints: "",
		bpipe_version: "tested with bpipe 0.9.8.7",
		author: "Sergi Sayols"

	output.dir = RNAtypes_OUT.replaceFirst("out=", "")
	def RNAtypes_FLAGS = RNAtypes_FOLDER   + " " +
	                     RNAtypes_PATTERN  + " " +
	                     RNAtypes_GTF      + " " +
	                     RNAtypes_OUT      + " " +
	                     RNAtypes_PRE      + " " +
	                     RNAtypes_SUF      + " " +
	                     RNAtypes_PAIRED   + " " +
	                     RNAtypes_STRANDED + " " +
	                     RNAtypes_MULTIMAP + " " +
	                     RNAtypes_FTYPE    + " " +
	                     RNAtypes_FTYPECOL + " " +
					     RNAtypes_CORES    + " " +
                         RNAtypes_EXTRA

	// run the chunk
	produce(RNAtypes_OUTNAME + ".counts.raw.png",
	        RNAtypes_OUTNAME + ".counts.per.png",
	        RNAtypes_OUTNAME + ".counts.rpk.png") {
		exec """
			export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES &&
			source ${TOOL_R}/env.sh &&
			if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi &&
			
			echo 'VERSION INFO'  1>&2 ;
			echo \$(${TOOL_R}/bin/Rscript --version 2>&1 | cut -d' ' -f5) 1>&2 ;
			echo '/VERSION INFO'  1>&2 ;
			
			${TOOL_R}/bin/Rscript ${TOOL_RNAtypes}/RNAtypes.R $RNAtypes_FLAGS &&
			mv ${RNAtypes_OUTNAME}.counts.* $output.dir
		""","RNAtypes"
	}
}

