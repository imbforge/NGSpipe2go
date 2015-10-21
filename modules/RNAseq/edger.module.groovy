//rule for task DE_edgeR from catalog RNAseq, version 1
//desc: Differential expression analysis using linears models and edgeR
DE_edgeR = {
	doc title: "DE_edgeR",
		desc:  "Differential expression analysis using linears models and edgeR",
		constraints: "",
		author: "Sergi Sayols"

	output.dir = DE_edgeR_OUTDIR
	def DE_edgeR_FLAGS = " targets="    + DE_edgeR_TARGETS + 
	                     " contrasts="  + DE_edgeR_CONTRASTS +
	                     " mmatrix="    + DE_edgeR_MMATRIX +
	                     " filter="     + DE_edgeR_FILTER +
	                     " prefix="     + DE_edgeR_PREFIX +
	                     " suffix="     + DE_edgeR_SUFFIX +
	                     " cwd="        + DE_edgeR_CWD +
	                     " robust="     + DE_edgeR_ROBUST +
                             " out="        + DE_edgeR_OUTDIR + "/DE_edgeR" +
                             " gtf="        + DE_edgeR_GTF

	// run the chunk
	produce("DE_edgeR.RData") {
		exec """
			export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES &&
			source ${TOOL_R}/env.sh &&
			if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi &&
			
			echo 'VERSION INFO'  1>&2 ;
			${TOOL_R}/bin/Rscript --version 1>&2 ;
			echo '/VERSION INFO' 1>&2 ;

			${TOOL_R}/bin/Rscript ${TOOL_EDGER}/DE_edgeR.R $DE_edgeR_FLAGS
		""","DE_edgeR"
	}
	
	//${TOOL_DEPENDENCIES}/R/default/bin/Rscript ${TOOL_DEPENDENCIES}/imb-forge/DE_edgeR/DE_edgeR.R $DE_edgeR_FLAGS
	
	forward input
}

