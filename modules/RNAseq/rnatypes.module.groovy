//rule for task RNAtypes from catalog RNAseq, version 1
//desc: Analysis of duplication rate on RNAseq analysis
RNAtypes = {
	doc title: "RNAtypes",
		desc:  "analysis of duplication rate on RNAseq analysis",
		constraints: "",
		bpipe_version: "tested with bpipe 0.9.8.7",
		author: "Sergi Sayols"

	output.dir = RNAtypes_OUTDIR.replaceFirst("out=", "")
	def RNAtypes_FLAGS = RNAtypes_FOLDER   + " " +
	                     RNAtypes_PATTERN  + " " +
	                     RNAtypes_GTF      + " " +
	                     RNAtypes_OUTDIR   + "/" + RNAtypes_OUTNAME + " " +
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
			
			echo 'VERSION INFO'  1>&2 ;
			echo \$(Rscript --version 2>&1 | cut -d' ' -f5) 1>&2 ;
			echo '/VERSION INFO'  1>&2 ;
			
			Rscript ${TOOL_RNAtypes}/RNAtypes.R $RNAtypes_FLAGS
		""","RNAtypes"
	}
}

