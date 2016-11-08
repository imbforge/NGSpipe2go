MappingStats = {
		doc title: "Statistics of mapping efficiency",
		desc:  "Counts the number of reads in the mapped bam, including total, unique, and mapped. Returns a plot of the results.",
		constraints: "Bam files produced by Bowtie 1.x. Might not work for other mappers/versions.",
		author: "Antonio Domingues, Anke Busch"

		output.dir = MAPPING_STATS_PLOTDIR

		produce(MAPPING_STATS_PLOTDIR + "/totalReads.pdf", MAPPING_STATS_PLOTDIR + "/totalReads.png") {

			exec """
				export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES &&
				source ${TOOL_R}/env.sh &&
				source ${TOOL_SAMTOOLS}/env.sh &&

			        echo 'VERSION INFO'  1>&2 &&
			        echo \$(Rscript --version 2>&1 | cut -d' ' -f5) 1>&2 &&
			        echo '/VERSION INFO'  1>&2 &&

			 	Rscript ${MAPPING_STATS_TOOL} ${MAPPING_STATS_DATADIR} ${MAPPING_STATS_PLOTDIR} ${ESSENTIAL_SAMPLE_PREFIX}

			""","MappingStats"
		}
}
