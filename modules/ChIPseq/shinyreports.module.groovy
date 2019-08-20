// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/ChIPseq/shinyreports.vars.groovy"

shinyReports = {
    doc title: "shinyReports",
        desc:  "creates the source code to compile the shiny and markdown reports",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols"

    output.dir = REPORTS

    def PREAMBLE = get_preamble("shinyReports")

    produce("shinyReports.txt") {
        exec """
            ${PREAMBLE} &&

            cp ${PIPELINE_ROOT}/tools/reports/shiny_chipseq_reporting_tool/server.R ${REPORTS}                &&
            cp ${PIPELINE_ROOT}/tools/reports/shiny_chipseq_reporting_tool/ui.R ${REPORTS}                    &&
            cp ${PIPELINE_ROOT}/tools/reports/shiny_chipseq_reporting_tool/ChIP.shinyrep.helpers.R ${REPORTS} &&
            cp ${PIPELINE_ROOT}/tools/reports/shiny_chipseq_reporting_tool/bustard.pl ${REPORTS}              &&
            cp ${PIPELINE_ROOT}/tools/reports/shiny_chipseq_reporting_tool/BustardSummary.toMD.xsl ${REPORTS} &&
            cp ${PIPELINE_ROOT}/tools/reports/shiny_chipseq_reporting_tool/styles.css ${REPORTS}              &&

            if [ -e "${REPORTS}/DEreport.Rmd" ]; then
                echo 'DEreport.Rmd already exists. Older copy will be kept and not overwritten';
            else
                cp ${PIPELINE_ROOT}/tools/reports/shiny_chipseq_reporting_tool/ChIPreport.Rmd ${REPORTS};
            fi &&

            PROJECT=\$(basename ${SHINYREPS_PROJECT})                              &&
            sed -i "2,2s/SHINYREPS_PROJECT/\${PROJECT}/" ${REPORTS}/ChIPreport.Rmd &&

            echo "SHINYREPS_PROJECT=${SHINYREPS_PROJECT}" >  $output &&
            echo "SHINYREPS_LOG=${SHINYREPS_LOG}"         >> $output &&
            echo "SHINYREPS_QC=${SHINYREPS_QC}"           >> $output &&
            echo "SHINYREPS_RES=${SHINYREPS_RES}"         >> $output &&
            echo "SHINYREPS_TARGETS=${SHINYREPS_TARGETS}" >> $output &&
            echo "SHINYREPS_BOWTIE_LOG=${SHINYREPS_BOWTIE_LOG}"      >> $output &&
            echo "SHINYREPS_BAMINDEX_LOG=${SHINYREPS_BAMINDEX_LOG}"  >> $output &&
            echo "SHINYREPS_MARKDUPS_LOG=${SHINYREPS_MARKDUPS_LOG}"  >> $output &&
            echo "SHINYREPS_EXTEND_LOG=${SHINYREPS_EXTEND_LOG}"      >> $output &&
            echo "SHINYREPS_FASTQC=${SHINYREPS_FASTQC}"   >> $output &&
            echo "SHINYREPS_FASTQC_LOG=${SHINYREPS_FASTQC_LOG}"   >> $output &&
            echo "SHINYREPS_IPSTRENGTH=${SHINYREPS_IPSTRENGTH}"      >> $output &&
            echo "SHINYREPS_IPSTRENGTH_LOG=${SHINYREPS_IPSTRENGTH_LOG}"      >> $output &&
            echo "SHINYREPS_PBC=${SHINYREPS_PBC}"         >> $output &&
            echo "SHINYREPS_PHANTOMPEAK=${SHINYREPS_PHANTOMPEAK}"    >> $output &&
            echo "SHINYREPS_PHANTOM_LOG=${SHINYREPS_PHANTOM_LOG}"    >> $output &&
            echo "SHINYREPS_BUSTARD=${SHINYREPS_BUSTARD}" >> $output &&
            echo "SHINYREPS_MACS2=${SHINYREPS_MACS2}"     >> $output &&
            echo "SHINYREPS_MACS2_LOG=${SHINYREPS_MACS2_LOG}"         >> $output &&
            echo "SHINYREPS_BLACKLIST_FILTER=${SHINYREPS_BLACKLIST_FILTER}" >> $output &&
            echo "SHINYREPS_PREFIX=${SHINYREPS_PREFIX}"   >> $output &&
            echo "SHINYREPS_PLOTS_COLUMN=${SHINYREPS_PLOTS_COLUMN}" >> $output &&
            echo "SHINYREPS_PEAK_ANNOTATION=${SHINYREPS_PEAK_ANNOTATION}" >> $output &&
            echo "SHINYREPS_GREAT=${SHINYREPS_GREAT}" >> $output &&
            echo "SHINYREPS_TRACKHUB_DONE=${SHINYREPS_TRACKHUB_DONE}" >> $output &&
            echo "SHINYREPS_TOOL_VERSIONS=${SHINYREPS_TOOL_VERSIONS}" >> $output
        ""","shinyReports"
    }
}

