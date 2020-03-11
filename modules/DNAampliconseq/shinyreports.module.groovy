// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/DNAampliconseq/shinyreports.vars.groovy"

shinyReports = {
    doc title: "shinyReports",
        desc:  "creates the source code to compile the shiny and markdown reports",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Sergi Sayols, modified by Frank Rühle"

    output.dir = REPORTS

    def PREAMBLE = get_preamble("shinyReports")

    produce("shinyReports.txt") {
        exec """
            ${PREAMBLE} &&

            cp ${PIPELINE_ROOT}/tools/reports/shiny_mps_reporting_tool/server.R ${REPORTS}                &&
            cp ${PIPELINE_ROOT}/tools/reports/shiny_mps_reporting_tool/ui.R ${REPORTS}                    &&
            cp ${PIPELINE_ROOT}/tools/reports/shiny_mps_reporting_tool/mps.shinyrep.helpers.R ${REPORTS}   &&
            cp ${PIPELINE_ROOT}/tools/reports/shiny_mps_reporting_tool/bustard.pl ${REPORTS}              &&
            cp ${PIPELINE_ROOT}/tools/reports/shiny_mps_reporting_tool/BustardSummary.toMD.xsl ${REPORTS} &&
            cp ${PIPELINE_ROOT}/tools/reports/shiny_mps_reporting_tool/styles.css ${REPORTS}              &&

            if [ -e "${REPORTS}/mps.report.Rmd" ]; then
                echo 'mps.report.Rmd already exists. Older copy will be kept and not overwritten';
            else
                cp ${PIPELINE_ROOT}/tools/reports/shiny_mps_reporting_tool/mps.report.Rmd ${REPORTS};
            fi &&

            PROJECT=\$(basename ${shinyReports_vars.project})                     &&
            sed -i "2,2s/SHINYREPS_PROJECT/\${PROJECT}/" ${REPORTS}/mps.report.Rmd &&

            echo "SHINYREPS_PROJECT=${shinyReports_vars.project}" >  $output &&
            echo "SHINYREPS_EXPDESIGN=${shinyReports_vars.expdesign}" >  $output &&
            echo "SHINYREPS_LOG=${shinyReports_vars.log}"         >> $output &&
            echo "SHINYREPS_QC=${shinyReports_vars.qc}"           >> $output &&
            echo "SHINYREPS_RES=${shinyReports_vars.res}"         >> $output &&
            echo "SHINYREPS_FASTQC_OUT=${shinyReports_vars.fastqc_out}"     >> $output &&
            echo "SHINYREPS_FASTQC_LOG=${shinyReports_vars.fastqc_log}"     >> $output &&

            echo "SHINYREPS_PEAR_OUT=${shinyReports_vars.pear_out}"     >> $output &&
            echo "SHINYREPS_PEAR_LOG=${shinyReports_vars.pear_log}"     >> $output &&
            echo "SHINYREPS_UMIEXTRACT_OUT=${shinyReports_vars.umiextract_out}"     >> $output &&
            echo "SHINYREPS_UMIEXTRACT_LOG=${shinyReports_vars.umiextract_log}"     >> $output &&
            echo "SHINYREPS_BARCODE_COUNT_OUT=${shinyReports_vars.barcode_count_out}"     >> $output &&
            echo "SHINYREPS_BARCODE_COUNT_LOG=${shinyReports_vars.barcode_count_log}"     >> $output &&

            echo "SHINYREPS_BAMINDEX_LOG=${shinyReports_vars.bamindex_log}" >> $output &&
            echo "SHINYREPS_BUSTARD=${shinyReports_vars.bustard}"         >> $output &&
            echo "SHINYREPS_BAM2BW_LOG=${shinyReports_vars.bam2bw_log}"         >> $output &&
            echo "SHINYREPS_PLOTS_COLUMN=${shinyReports_vars.plots_column}"     >> $output &&
            echo "SHINYREPS_INSERTSIZE=${shinyReports_vars.insertsize}" >> $output &&
            echo "SHINYREPS_TRACKHUB_DONE=${shinyReports_vars.trackhub_done}" >> $output &&
            echo "SHINYREPS_TOOL_VERSIONS=${shinyReports_vars.tool_versions}" >> $output &&
            echo "SHINYREPS_TARGET=${shinyReports_vars.target}"   >> $output 
        ""","shinyReports"
    }
}
