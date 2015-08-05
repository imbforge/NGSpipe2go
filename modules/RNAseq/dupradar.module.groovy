dupRadar = {
	doc title: "dupRadar",
		desc:  "analysis of duplication rate on RNAseq analysis",
		constraints: "",
		author: "Sergi Sayols"

	output.dir = DUPRADAR_OUTDIR
	def DUPRADAR_FLAGS = " gtf="      + DUPRADAR_GTF      +
	                     " stranded=" + DUPRADAR_STRANDED + 
			    		 " paired="   + DUPRADAR_PAIRED   +
						 " outdir="   + DUPRADAR_OUTDIR   +
				    	 " threads="  + DUPRADAR_THREADS

	// run the chunk
	transform(".bam") to("_dupRadar.png") {
		exec """
			export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES &&
			source ${TOOL_R}/env.sh &&
			if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi &&

			${TOOL_R}/bin/Rscript ${TOOL_DUPRADAR}/dupRadar.R $input $DUPRADAR_FLAGS
		""","dupRadar"
	}

	forward input
}

