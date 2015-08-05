//rule for task RNAtypes from catalog RNAseq, version 1
//desc: Analysis of duplication rate on RNAseq analysis
RNAtypes = {
	doc title: "RNAtypes",
		desc:  "analysis of duplication rate on RNAseq analysis",
		constraints: "",
		author: "Sergi Sayols"

	output.dir = QC + "/RNAtypes/"
	def RNAtypes_FLAGS = " folder="   + RNAtypes_FOLDER   +
	                     " pattern="  + RNAtypes_PATTERN  +
	                     " gtf="      + RNAtypes_GTF +
	                     " out="      + RNAtypes_OUT +
	                     " pre="      + RNAtypes_PRE +
	                     " suf="      + RNAtypes_SUF +
	                     " paired="   + RNAtypes_PAIRED   +
	                     " stranded=" + RNAtypes_STRANDED +
	                     " multimap=" + RNAtypes_MULTIMAP +
	                     " ftype="    + RNAtypes_FTYPE    +
	                     " ftypecol=" + RNAtypes_FTYPECOL +
					     " cores="    + RNAtypes_CORES

	// run the chunk
	produce(RNAtypes_OUT + ".counts.raw.png",
	        RNAtypes_OUT + ".counts.per.png",
	        RNAtypes_OUT + ".counts.rpk.png") {
		exec """
			export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES &&
			source ${TOOL_R}/env.sh &&
			if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi &&

			${TOOL_R}/bin/Rscript ${TOOL_RNAtypes}/RNAtypes.R $RNAtypes_FLAGS &&
			mv ${RNAtypes_OUT}.counts.* $output.dir
		""","RNAtypes"
	}
}

