//rule for task shinyReports from catalog ChIPseq, version 1
//desc: Creates the source code to compile the shiny and markdown reports
shinyReports = {
	doc title: "shinyReports",
		desc:  "creates the source code to compile the shiny and markdown reports",
		constraints: "",
		author: "Sergi Sayols"

	output.dir = REPORTS
	
	produce("shinyReports.txt") {
		exec """
			
			cp ${MODULE_FOLDER}/../tools/reports/shiny_chipseq_reporting_tool/server.R ${REPORTS}                &&
			cp ${MODULE_FOLDER}/../tools/reports/shiny_chipseq_reporting_tool/ui.R ${REPORTS}                    &&
			cp ${MODULE_FOLDER}/../tools/reports/shiny_chipseq_reporting_tool/ChIP.shinyrep.helpers.R ${REPORTS} &&
			cp ${MODULE_FOLDER}/../tools/reports/shiny_chipseq_reporting_tool/bustard.pl ${REPORTS}              &&
			cp ${MODULE_FOLDER}/../tools/reports/shiny_chipseq_reporting_tool/BustardSummary.toMD.xsl ${REPORTS} &&
			cp ${MODULE_FOLDER}/../tools/reports/shiny_rnaseq_reporting_tool/styles.css ${REPORTS}               &&
			
			if [ -e "${REPORTS}/DEreport.Rmd" ]; then
				echo 'DEreport.Rmd already exists. Older copy will be kept and not overwritten';
			else
				cp ${MODULE_FOLDER}/../tools/reports/shiny_chipseq_reporting_tool/ChIPreport.Rmd ${REPORTS};
			fi &&
			
	
			PROJECT=\$(basename ${SHINYREPS_PROJECT})                              &&
			sed -i "2,2s/SHINYREPS_PROJECT/\${PROJECT}/" ${REPORTS}/ChIPreport.Rmd &&
	
			echo "SHINYREPS_PROJECT=${SHINYREPS_PROJECT}" >  $output &&
			echo "SHINYREPS_LOG=${SHINYREPS_LOG}"         >> $output &&
			echo "SHINYREPS_QC=${SHINYREPS_QC}"           >> $output &&
			echo "SHINYREPS_RES=${SHINYREPS_RES}"         >> $output &&
			echo "SHINYREPS_TARGETS=${SHINYREPS_TARGETS}" >> $output &&
			echo "SHINYREPS_BOWTIE_LOG=${SHINYREPS_BOWTIE_LOG}"     >> $output &&
			echo "SHINYREPS_MARKDUPS_LOG=${SHINYREPS_MARKDUPS_LOG}" >> $output &&
			echo "SHINYREPS_FASTQC=${SHINYREPS_FASTQC}"   >> $output &&
			echo "SHINYREPS_IPSTRENGTH=${SHINYREPS_IPSTRENGTH}"     >> $output &&
			echo "SHINYREPS_PBC=${SHINYREPS_PBC}"         >> $output &&
			echo "SHINYREPS_PHANTOMPEAK=${SHINYREPS_PHANTOMPEAK}"   >> $output &&
			echo "SHINYREPS_BUSTARD=${SHINYREPS_BUSTARD}" >> $output &&
			echo "SHINYREPS_MACS2=${SHINYREPS_MACS2}"     >> $output &&
			echo "SHINYREPS_PREFIX=${SHINYREPS_PREFIX}"   >> $output &&
			echo "SHINYREPS_PLOTS_COLUMN=${SHINYREPS_PLOTS_COLUMN}" >> $output
		""","shinyReports"
	}
}

