CombinedStats = {
	doc title: "Summary of all trimming and filtering statistics of raw reads",
	desc:  "Counts the number of reads in the original reads file, after adapter trimming, after quality filtering and after duplicate removal and plots results",
	author: "Anke Busch"

	output.dir = COMBINED_STATS_PLOTDIR

	produce(COMBINED_STATS_PLOTDIR + "/allTrimmingStats.pdf", COMBINED_STATS_PLOTDIR + "/allTrimmingStats.png") {

		exec """
            module load R/${R_VERSION} &&
			
            Rscript ${COMBINED_STATS_PLOT_TOOL} ${COMBINED_STATS_DATADIR} ${COMBINED_STATS_PLOTDIR} ${ESSENTIAL_SAMPLE_PREFIX}

		""","CombinedStats"
	}
}
