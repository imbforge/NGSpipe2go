DedupStats = {
	doc title: "Statistics of unique reads",
		desc:  "Counts the number of reads in the original reads file, and after PCR duplicate removal, and plots results",
		author: "Antonio Domingues"

   PLOT_TOOL = PLOT_TOOL_PATH + "/PCRDuplicatesPlot.R"

	produce(REMOVE_DUP_OUTDIR + "/figure/PCRDuplicates.pdf",
           REMOVE_DUP_OUTDIR + "/figure/PCRDuplicates.png") {

      exec """
         if [ -n "\$LSB_JOBID" ]; then
            export TMPDIR=/jobdir/\${LSB_JOBID};
         fi &&

         echo 'VERSION INFO'  1>&2 &&
         echo \$(Rscript --version 2>&1 | cut -d' ' -f5) 1>&2 &&
         echo '/VERSION INFO'  1>&2 &&

         cd ${REMOVE_DUP_OUTDIR} &&
			for i in *.fastq2table; do wc -l $i >> dedup.stats.txt; done &&
         for i in *.unique; do wc -l $i >> dedup.stats.txt; done &&
         Rscript ${PLOT_TOOL}

		""","DedupStats"
	}
}
