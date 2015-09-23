//rule for task shinyReports from catalog miscellaneous, version 1
//desc: creates the source code to compile the shiny and markdown reports
shinyReports = {
	doc title: "shinyReports",
		desc:  "creates the source code to compile the shiny and markdown reports",
		constraints: "",
		author: "Sergi Sayols"

	exec """
		if [ ! -d "${REPORTS}" ]; then
			mkdir -p ${REPORTS};
		fi &&
		
		cp ${TOOL_DEPENDENCIES}/imb-forge/reports/rnaseq/server.R ${REPORTS}                &&
		cp ${TOOL_DEPENDENCIES}/imb-forge/reports/rnaseq/ui.R ${REPORTS}                    &&
		cp ${TOOL_DEPENDENCIES}/imb-forge/reports/rnaseq/DE.shinyrep.helpers.R ${REPORTS}   &&
		cp ${TOOL_DEPENDENCIES}/imb-forge/reports/rnaseq/bustard.pl ${REPORTS}              &&
		cp ${TOOL_DEPENDENCIES}/imb-forge/reports/rnaseq/BustardSummary.toMD.xsl ${REPORTS} &&
		
		if [ -e "${REPORTS}/DEreport.Rmd" ]; then
			echo 'DEreport.Rmd already exists. Older copy will be kept and not overwritten';
		else
			cp ${TOOL_DEPENDENCIES}/imb-forge/reports/rnaseq/DEreport.Rmd ${REPORTS};
		fi &&
		
		PROJECT=\$(basename ${SHINYREPS_PROJECT})                            &&
		sed -i "2,2s/SHINYREPS_PROJECT/\${PROJECT}/" ${REPORTS}/DEreport.Rmd &&
		
		echo "SHINYREPS_PROJECT=${SHINYREPS_PROJECT}" >  ${REPORTS}/shinyReports.txt &&
		echo "SHINYREPS_ORG=${SHINYREPS_ORG}"         >> ${REPORTS}/shinyReports.txt &&
		echo "SHINYREPS_DB=${SHINYREPS_DB}"           >> ${REPORTS}/shinyReports.txt &&
		echo "SHINYREPS_LOG=${SHINYREPS_LOG}"         >> ${REPORTS}/shinyReports.txt &&
		echo "SHINYREPS_QC=${SHINYREPS_QC}"           >> ${REPORTS}/shinyReports.txt &&
		echo "SHINYREPS_RES=${SHINYREPS_RES}"         >> ${REPORTS}/shinyReports.txt &&
		echo "SHINYREPS_PREFIX=${SHINYREPS_PREFIX}"   >> ${REPORTS}/shinyReports.txt &&
		echo "SHINYREPS_STAR_LOG=${SHINYREPS_STAR_LOG}"       >> ${REPORTS}/shinyReports.txt &&
		echo "SHINYREPS_STAR_SUFFIX=${SHINYREPS_STAR_SUFFIX}" >> ${REPORTS}/shinyReports.txt &&
		echo "SHINYREPS_STARparms_SUFFIX=${SHINYREPS_STARparms_SUFFIX}" >> ${REPORTS}/shinyReports.txt &&
		echo "SHINYREPS_FASTQC_LOG=${SHINYREPS_FASTQC_LOG}"   >> ${REPORTS}/shinyReports.txt &&
		echo "SHINYREPS_DUPRADAR_LOG=${SHINYREPS_DUPRADAR_LOG}" >> ${REPORTS}/shinyReports.txt &&
		echo "SHINYREPS_RNATYPES_LOG=${SHINYREPS_RNATYPES_LOG}" >> ${REPORTS}/shinyReports.txt &&
		echo "SHINYREPS_GENEBODYCOV_LOG=${SHINYREPS_GENEBODYCOV_LOG}" >> ${REPORTS}/shinyReports.txt &&
		echo "SHINYREPS_BUSTARD=${SHINYREPS_BUSTARD}"         >> ${REPORTS}/shinyReports.txt &&
		echo "SHINYREPS_DE_EDGER=${SHINYREPS_DE_EDGER}"       >> ${REPORTS}/shinyReports.txt &&
		echo "SHINYREPS_SUBREAD=${SHINYREPS_SUBREAD}"       >> ${REPORTS}/shinyReports.txt &&
		echo "SHINYREPS_SUBREAD_SUFFIX=${SHINYREPS_SUBREAD_SUFFIX}"       >> ${REPORTS}/shinyReports.txt
	""","shinyReports"
}

