DE_DESeq2 = {
	doc title: "DE_DESeq2",
		desc:  "Differential expression analysis using linears models and DESeq2",
		constraints: "Only simple contrasts in 1-factor design. Include always the intercept. Always gene filtering",
		author: "Sergi Sayols"

	output.dir = DE_DESeq2_OUTDIR
	def DE_DESeq2_FLAGS = " targets="  + DE_DESeq2_TARGETS + 
	                     " contrasts=" + DE_DESeq2_CONTRASTS +
	                     " mmatrix="   + DE_DESeq2_MMATRIX +
	                     " filter="    + DE_DESeq2_FILTER +
	                     " prefix="    + DE_DESeq2_PREFIX +
	                     " suffix="    + DE_DESeq2_SUFFIX +
	                     " cwd="       + DE_DESeq2_CWD +
                         " out="       + DE_DESeq2_OUTDIR + "/DE_DESeq2"

	// run the chunk
	produce("DE_DESeq2.RData") {
		exec """
			export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES &&
			source ${TOOL_R}/env.sh &&
			if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi &&
			
			echo '<<<<<<<<<<<<<<<<' &&
			echo ${TOOL_R}/bin/Rscript ${TOOL_DESeq2}/DE_DESeq2.R $DE_DESeq2_FLAGS &&
			echo '>>>>>>>>>>>>>>>>' &&
			
			${TOOL_R}/bin/Rscript ${TOOL_DESeq2}/DE_DESeq2.R $DE_DESeq2_FLAGS
		""","DE_DESeq2"
	}

	forward input
}

