DedupStats = {
	doc title: "Statistics of unique reads",
		desc:  "Counts the number of reads in the original reads file, and after PCR duplicate removal, and plots results",
		author: "Antonio Domingues"


	produce(REMOVE_DUP_OUTDIR + "/figure/PCRDuplicates.pdf",
           REMOVE_DUP_OUTDIR + "/figure/PCRDuplicates.png") {

      exec """
         echo 'VERSION INFO'  1>&2 &&
         echo \$(${TOOL_R}/bin/Rscript --version 2>&1 | cut -d' ' -f5) 1>&2 &&
         echo '/VERSION INFO'  1>&2 &&

         cd ${REMOVE_DUP_OUTDIR} &&
         ${TOOL_R}/bin/Rscript ${DEDUP_PLOT_TOOL}

		""","DedupStats"
	}
}
