diffbind = {
	doc title: "diffbind",
		desc:  "Differential binding analysis using Diffbind",
		constraints: "",
		bpipe_version: "tested with bpipe 0.9.8.7",
		author: "Sergi Sayols"

	output.dir = RESULTS + "/diffbind"
	def DIFFBIND_FLAGS = DIFFBIND_TARGETS   + " " + 
                         DIFFBIND_CONTRASTS + " " +
                         DIFFBIND_CWD       + " " +
                         DIFFBIND_OUTDIR    + " " +
                         DIFFBIND_BAMS      + " " +
                         DIFFBIND_FRAGSIZE  + " " +
                         DIFFBIND_ANNOTATE  + " " +
                         DIFFBIND_TSS       + " " +
                         DIFFBIND_TXDB      + " " +
                         DIFFBIND_ANNODB    + " " +
                         DIFFBIND_EXTRA

	// run the chunk
	produce("diffbind.pdf", "diffbind.xls") {
		exec """
			module load R/${R_VERSION} &&
			
			Rscript ${TOOL_DIFFBIND}/diffbind.R $DIFFBIND_FLAGS
		""","diffbind"
	}

	forward input
}

