//rule for task DE_edgeR from catalog RNAseq, version 1
//desc: Differential expression analysis using linears models and edgeR
DE_edgeR = {
	doc title: "DE_edgeR",
		desc:  "Differential expression analysis using linears models and edgeR",
		constraints: "",
		bpipe_version: "tested with bpipe 0.9.8.7",
		author: "Sergi Sayols"

	output.dir = DE_edgeR_OUTDIR.replcaeFirst("out=", "")
	def DE_edgeR_FLAGS = DE_edgeR_TARGETS   + " " + 
	                     DE_edgeR_CONTRASTS + " " +
	                     DE_edgeR_MMATRIX   + " " +
	                     DE_edgeR_FILTER    + " " +
	                     DE_edgeR_PREFIX    + " " +
	                     DE_edgeR_SUFFIX    + " " +
	                     DE_edgeR_CWD       + " " +
	                     DE_edgeR_ROBUST    + " " +
                         DE_edgeR_GTF       + " " +
                         DE_edgeR_OUTDIR    + "/DE_edgeR " +
                         DE_edgeR_EXTRA

	// run the chunk
	produce("DE_edgeR.RData") {
		exec """
			module load R/${R_VERSION} &&
			
		    Rscript ${TOOL_EDGER}/DE_edgeR.R $DE_edgeR_FLAGS
		""","DE_edgeR"
	}
	
	forward input
}

