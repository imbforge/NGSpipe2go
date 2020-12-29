DedupStats = {
	doc title: "Statistics of unique reads",
		desc:  "Counts the number of reads in the original reads file, and after PCR duplicate removal, and plots results",
		author: "Antonio Domingues"


	produce(REMOVE_DUP_OUTDIR + "/figure/PCRDuplicates.pdf",
           REMOVE_DUP_OUTDIR + "/figure/PCRDuplicates.png") {

      exec """
         module load R/${R_VERSION} &&

         cd ${REMOVE_DUP_OUTDIR} &&
         Rscript ${DEDUP_PLOT_TOOL}

		""","DedupStats"
	}
}
