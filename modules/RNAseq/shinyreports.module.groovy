//rule for task shinyReports from catalog miscellaneous, version 1
//desc: creates the source code to compile the shiny and markdown reports
shinyReports = {
	doc title: "shinyReports",
		desc:  "creates the source code to compile the shiny and markdown reports",
		constraints: "",
		bpipe_version: "tested with bpipe 0.9.8.7",
		author: "Sergi Sayols"
	
	output.dir = REPORTS
	
	produce("shinyReports.txt") {
		exec """
			
			cp ${MODULE_FOLDER}/../tools/reports/shiny_rnaseq_reporting_tool/server.R ${REPORTS}                &&
			cp ${MODULE_FOLDER}/../tools/reports/shiny_rnaseq_reporting_tool/ui.R ${REPORTS}                    &&
			cp ${MODULE_FOLDER}/../tools/reports/shiny_rnaseq_reporting_tool/DE.shinyrep.helpers.R ${REPORTS}   &&
			cp ${MODULE_FOLDER}/../tools/reports/shiny_rnaseq_reporting_tool/bustard.pl ${REPORTS}              &&
			cp ${MODULE_FOLDER}/../tools/reports/shiny_rnaseq_reporting_tool/BustardSummary.toMD.xsl ${REPORTS} &&
			cp ${MODULE_FOLDER}/../tools/reports/shiny_rnaseq_reporting_tool/styles.css ${REPORTS}              &&
			
			if [ -e "${REPORTS}/DEreport.Rmd" ]; then
				echo 'DEreport.Rmd already exists. Older copy will be kept and not overwritten';
			else
				cp ${MODULE_FOLDER}/../tools/reports/shiny_rnaseq_reporting_tool/DEreport.Rmd ${REPORTS};
			fi &&
			
			PROJECT=\$(basename ${SHINYREPS_PROJECT})                            &&
			sed -i "2,2s/SHINYREPS_PROJECT/\${PROJECT}/" ${REPORTS}/DEreport.Rmd &&
			
			echo "SHINYREPS_PROJECT=${SHINYREPS_PROJECT}" >  $output &&
			echo "SHINYREPS_ORG=${SHINYREPS_ORG}"         >> $output &&
			echo "SHINYREPS_DB=${SHINYREPS_DB}"           >> $output &&
			echo "SHINYREPS_LOG=${SHINYREPS_LOG}"         >> $output &&
			echo "SHINYREPS_QC=${SHINYREPS_QC}"           >> $output &&
			echo "SHINYREPS_RES=${SHINYREPS_RES}"         >> $output &&
			echo "SHINYREPS_PREFIX=${SHINYREPS_PREFIX}"   >> $output &&
			echo "SHINYREPS_STAR_LOG=${SHINYREPS_STAR_LOG}"       >> $output &&
			echo "SHINYREPS_STAR_SUFFIX=${SHINYREPS_STAR_SUFFIX}" >> $output &&
			echo "SHINYREPS_STARparms_SUFFIX=${SHINYREPS_STARparms_SUFFIX}" >> $output &&
			echo "SHINYREPS_FASTQC_OUT=${SHINYREPS_FASTQC_OUT}"     >> $output &&
			echo "SHINYREPS_FASTQC_LOG=${SHINYREPS_FASTQC_LOG}"     >> $output &&
			echo "SHINYREPS_BAMINDEX_LOG=${SHINYREPS_BAMINDEX_LOG}"     >> $output &&
			echo "SHINYREPS_DUPRADAR_LOG=${SHINYREPS_DUPRADAR_LOG}" >> $output &&
			echo "SHINYREPS_RNATYPES_LOG=${SHINYREPS_RNATYPES_LOG}" >> $output &&
			echo "SHINYREPS_GENEBODYCOV_LOG=${SHINYREPS_GENEBODYCOV_LOG}" >> $output &&
			echo "SHINYREPS_BUSTARD=${SHINYREPS_BUSTARD}"         >> $output &&
			echo "SHINYREPS_DE_EDGER=${SHINYREPS_DE_EDGER}"       >> $output &&
			echo "SHINYREPS_DE_DESEQ=${SHINYREPS_DE_DESEQ}"       >> $output &&
			echo "SHINYREPS_DE_DESEQ_MM=${SHINYREPS_DE_DESEQ_MM}" >> $output &&
			echo "SHINYREPS_SUBREAD=${SHINYREPS_SUBREAD}"         >> $output &&
			echo "SHINYREPS_SUBREAD_SUFFIX=${SHINYREPS_SUBREAD_SUFFIX}"  >> $output &&
			echo "SHINYREPS_SUBREAD_LOG=${SHINYREPS_SUBREAD_LOG}"        >> $output &&
			echo "SHINYREPS_BAM2BW_LOG=${SHINYREPS_BAM2BW_LOG}"          >> $output &&
			echo "SHINYREPS_MARKDUPS_LOG=${SHINYREPS_MARKDUPS_LOG}"      >> $output &&
			echo "SHINYREPS_EDGER_LOGS=${SHINYREPS_EDGER_LOGS}"          >> $output &&
			echo "SHINYREPS_DESEQ_LOGS=${SHINYREPS_DESEQ_LOGS}"          >> $output &&
			echo "SHINYREPS_PLOTS_COLUMN=${SHINYREPS_PLOTS_COLUMN}"      >> $output &&
			echo "SHINYREPS_INFEREXPERIMENT_LOGS=${SHINYREPS_INFEREXPERIMENT_LOGS}" >> $output &&
			echo "SHINYREPS_QUALIMAP_LOGS=${SHINYREPS_QUALIMAP_LOGS}" >> $output
		""","shinyReports"
	}
}

