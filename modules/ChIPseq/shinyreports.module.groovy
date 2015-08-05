//rule for task shinyReports from catalog ChIPseq, version 1
//desc: Creates the source code to compile the shiny and markdown reports
shinyReports = {
	doc title: "shinyReports",
		desc:  "creates the source code to compile the shiny and markdown reports",
		constraints: "",
		author: "Sergi Sayols"

	exec """
		if [ ! -d "${REPORTS}" ]; then
			mkdir -p ${REPORTS};
		fi &&

		cp ${TOOL_DEPENDENCIES}/imb-forge/reports/chipseq/server.R ${REPORTS}                &&
		cp ${TOOL_DEPENDENCIES}/imb-forge/reports/chipseq/ui.R ${REPORTS}                    &&
		cp ${TOOL_DEPENDENCIES}/imb-forge/reports/chipseq/ChIP.shinyrep.helpers.R ${REPORTS} &&
		cp ${TOOL_DEPENDENCIES}/imb-forge/reports/chipseq/bustard.pl ${REPORTS}              &&
		cp ${TOOL_DEPENDENCIES}/imb-forge/reports/chipseq/BustardSummary.toMD.xsl ${REPORTS} &&
		cp ${TOOL_DEPENDENCIES}/imb-forge/reports/chipseq/ChIPreport.Rmd ${REPORTS}          &&

		PROJECT=\$(basename ${SHINYREPS_PROJECT})                              &&
		sed -i "2,2s/SHINYREPS_PROJECT/\${PROJECT}/" ${REPORTS}/ChIPreport.Rmd &&

		echo "SHINYREPS_PROJECT=${SHINYREPS_PROJECT}" >  ${REPORTS}/shinyReports.txt &&
		echo "SHINYREPS_LOG=${SHINYREPS_LOG}"         >> ${REPORTS}/shinyReports.txt &&
		echo "SHINYREPS_QC=${SHINYREPS_QC}"           >> ${REPORTS}/shinyReports.txt &&
		echo "SHINYREPS_RES=${SHINYREPS_RES}"         >> ${REPORTS}/shinyReports.txt &&
		echo "SHINYREPS_TARGETS=${SHINYREPS_TARGETS}" >> ${REPORTS}/shinyReports.txt &&
		echo "SHINYREPS_BOWTIE_LOG=${SHINYREPS_BOWTIE_LOG}"     >> ${REPORTS}/shinyReports.txt &&
		echo "SHINYREPS_MARKDUPS_LOG=${SHINYREPS_MARKDUPS_LOG}" >> ${REPORTS}/shinyReports.txt &&
		echo "SHINYREPS_FASTQC=${SHINYREPS_FASTQC}"   >> ${REPORTS}/shinyReports.txt &&
		echo "SHINYREPS_IPSTRENGTH=${SHINYREPS_IPSTRENGTH}"     >> ${REPORTS}/shinyReports.txt &&
		echo "SHINYREPS_PBC=${SHINYREPS_PBC}"         >> ${REPORTS}/shinyReports.txt &&
		echo "SHINYREPS_PHANTOMPEAK=${SHINYREPS_PHANTOMPEAK}"   >> ${REPORTS}/shinyReports.txt &&
		echo "SHINYREPS_BUSTARD=${SHINYREPS_BUSTARD}" >> ${REPORTS}/shinyReports.txt &&
		echo "SHINYREPS_MACS2=${SHINYREPS_MACS2}"     >> ${REPORTS}/shinyReports.txt &&
		echo "SHINYREPS_PREFIX=${SHINYREPS_PREFIX}"   >> ${REPORTS}/shinyReports.txt
	""","shinyReports"
}

