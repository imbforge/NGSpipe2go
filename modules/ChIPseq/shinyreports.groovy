shinyReports = {
    doc title: "shinyReports",
        desc:  "creates the source code to compile the shiny and markdown reports",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
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

            PROJECT=\$(basename ${shinyReports_vars.project})                              &&
            sed -i "2,2s/SHINYREPS_PROJECT/\${PROJECT}/" ${REPORTS}/ChIPreport.Rmd &&

            echo "SHINYREPS_PROJECT=${shinyReports_vars.project}" >  $output &&
            echo "SHINYREPS_LOG=${shinyReports_vars.log}"         >> $output &&
            echo "SHINYREPS_QC=${shinyReports_vars.qc}"           >> $output &&
            echo "SHINYREPS_RES=${shinyReports_vars.res}"         >> $output &&
            echo "SHINYREPS_TARGETS=${shinyReports_vars.targets}" >> $output &&
            echo "SHINYREPS_PAIRED=${shinyReports_vars.paired}"      >> $output &&
            echo "SHINYREPS_BOWTIE_LOG=${shinyReports_vars.bowtie_log}"      >> $output &&
            echo "SHINYREPS_BAMINDEX_LOG=${shinyReports_vars.bamindex_log}"  >> $output &&
            echo "SHINYREPS_MARKDUPS_LOG=${shinyReports_vars.markdups_log}"  >> $output &&
            echo "SHINYREPS_EXTEND_LOG=${shinyReports_vars.extend_log}"      >> $output &&
            echo "SHINYREPS_FASTQC=${shinyReports_vars.fastqc}"   >> $output &&
            echo "SHINYREPS_FASTQC_LOG=${shinyReports_vars.fastqc_log}"   >> $output &&
            echo "SHINYREPS_INSERTSIZE=${shinyReports_vars.insertsize}" >> $output &&
            echo "SHINYREPS_IPSTRENGTH=${shinyReports_vars.ipstrength}"      >> $output &&
            echo "SHINYREPS_IPSTRENGTH_LOG=${shinyReports_vars.ipstrength_log}"      >> $output &&
            echo "SHINYREPS_PBC=${shinyReports_vars.pbc}"         >> $output &&
            echo "SHINYREPS_PHANTOMPEAK=${shinyReports_vars.phantompeak}"    >> $output &&
            echo "SHINYREPS_PHANTOM_LOG=${shinyReports_vars.phantom_log}"    >> $output &&
            echo "SHINYREPS_BUSTARD=${shinyReports_vars.bustard}" >> $output &&
            echo "SHINYREPS_MACS2=${shinyReports_vars.macs2}"     >> $output &&
            echo "SHINYREPS_MACS2_LOG=${shinyReports_vars.macs2_log}"         >> $output &&
            echo "SHINYREPS_BLACKLIST_FILTER=${shinyReports_vars.blacklist_filter}" >> $output &&
            echo "SHINYREPS_PREFIX=${shinyReports_vars.prefix}"   >> $output &&
            echo "SHINYREPS_PLOTS_COLUMN=${shinyReports_vars.plots_column}" >> $output &&
            echo "SHINYREPS_PEAK_ANNOTATION=${shinyReports_vars.peak_annotation}" >> $output &&
            echo "SHINYREPS_GREAT=${shinyReports_vars.great}" >> $output &&
            echo "SHINYREPS_DIFFBIND=${shinyReports_vars.diffbind}" >> $output           &&
            echo "SHINYREPS_TRACKHUB_DONE=${shinyReports_vars.trackhub_done}" >> $output &&
            echo "SHINYREPS_TOOL_VERSIONS=${shinyReports_vars.tool_versions}" >> $output
        ""","shinyReports"
    }
}

