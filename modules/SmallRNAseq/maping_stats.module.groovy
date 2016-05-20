MappingStatsPlot = {
	doc title: "Statistics of mapping efficiency",
		desc:  "Counts the number of reads in the mapped bam, including total, unique, and mapped. Returns a plot of the results.",
      constraints: "Bam files produced by Bowtie 1.x. Might not work for other mappers/versions.",
		author: "Antonio Domingues"

   MAPPING_STATS_TOOL = MAPPING_STATS_TOOL_PATH + "/mapping_stats_bowtie1.R"

	produce(MAPPING_STATS_OUTDIR + "/figure/totalReads.pdf") {

      exec """
         if [ -n "\$LSB_JOBID" ]; then
            export TMPDIR=/jobdir/\${LSB_JOBID};
         fi &&

         echo 'VERSION INFO'  1>&2 &&
         echo \$(${TOOL_R}/bin/Rscript --version 2>&1 | cut -d' ' -f5) 1>&2 &&
         echo '/VERSION INFO'  1>&2 &&

         cd ${MAPPING_STATS_OUTDIR} &&
			${TOOL_R}/bin/Rscript ${MAPPING_STATS_TOOL}

		""","MappingStatsPlot"
	}
}
