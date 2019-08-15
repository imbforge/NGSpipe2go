// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/scRNAseq/shinyreports.vars.groovy"

shinyReports = {
    doc title: "shinyReports",
        desc:  "creates the source code to compile the shiny and markdown reports",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols, modified by Frank Rühle"

    output.dir = REPORTS

    produce("shinyReports.txt") {
        exec """
            cp ${PIPELINE_ROOT}/tools/reports/shiny_scrnaseq_reporting_tool/server.R ${REPORTS}                &&
            cp ${PIPELINE_ROOT}/tools/reports/shiny_scrnaseq_reporting_tool/ui.R ${REPORTS}                    &&
            cp ${PIPELINE_ROOT}/tools/reports/shiny_scrnaseq_reporting_tool/sc.shinyrep.helpers.R ${REPORTS}   &&
            cp ${PIPELINE_ROOT}/tools/reports/shiny_scrnaseq_reporting_tool/bustard.pl ${REPORTS}              &&
            cp ${PIPELINE_ROOT}/tools/reports/shiny_scrnaseq_reporting_tool/BustardSummary.toMD.xsl ${REPORTS} &&
            cp ${PIPELINE_ROOT}/tools/reports/shiny_scrnaseq_reporting_tool/styles.css ${REPORTS}              &&

            if [ -e "${REPORTS}/sc.report.Rmd" ]; then
                echo 'sc.report.Rmd already exists. Older copy will be kept and not overwritten';
            else
                cp ${PIPELINE_ROOT}/tools/reports/shiny_scrnaseq_reporting_tool/sc.report.Rmd ${REPORTS};
            fi &&

            PROJECT=\$(basename ${SHINYREPS_PROJECT})                            &&
            sed -i "2,2s/SHINYREPS_PROJECT/\${PROJECT}/" ${REPORTS}/sc.report.Rmd &&

            echo "SHINYREPS_PROJECT=${SHINYREPS_PROJECT}" >  $output &&
            echo "SHINYREPS_ORG=${SHINYREPS_ORG}"         >> $output &&
            echo "SHINYREPS_DB=${SHINYREPS_DB}"           >> $output &&
            echo "SHINYREPS_LOG=${SHINYREPS_LOG}"         >> $output &&
            echo "SHINYREPS_PAIRED=${SHINYREPS_PAIRED}"   >> $output &&
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
            echo "SHINYREPS_RNATYPES=${SHINYREPS_RNATYPES}" >> $output &&
            echo "SHINYREPS_RNATYPES_SUFFIX=${SHINYREPS_RNATYPES_SUFFIX}" >> $output &&
            echo "SHINYREPS_GENEBODYCOV_LOG=${SHINYREPS_GENEBODYCOV_LOG}" >> $output &&
            echo "SHINYREPS_BUSTARD=${SHINYREPS_BUSTARD}"         >> $output &&
            echo "SHINYREPS_SUBREAD=${SHINYREPS_SUBREAD}"         >> $output &&
            echo "SHINYREPS_SUBREAD_SUFFIX=${SHINYREPS_SUBREAD_SUFFIX}"  >> $output &&
            echo "SHINYREPS_SUBREAD_LOG=${SHINYREPS_SUBREAD_LOG}"        >> $output &&
            echo "SHINYREPS_UMICOUNT=${SHINYREPS_UMICOUNT}"         >> $output &&
            echo "SHINYREPS_UMICOUNT_LOG=${SHINYREPS_UMICOUNT_LOG}"        >> $output &&
            echo "SHINYREPS_BAM2BW_LOG=${SHINYREPS_BAM2BW_LOG}"          >> $output &&
            echo "SHINYREPS_MARKDUPS_LOG=${SHINYREPS_MARKDUPS_LOG}"      >> $output &&
            echo "SHINYREPS_PLOTS_COLUMN=${SHINYREPS_PLOTS_COLUMN}"      >> $output &&
            echo "SHINYREPS_INFEREXPERIMENT_LOGS=${SHINYREPS_INFEREXPERIMENT_LOGS}" >> $output &&
            echo "SHINYREPS_QUALIMAP_LOGS=${SHINYREPS_QUALIMAP_LOGS}" >> $output &&
            echo "SHINYREPS_INSERTSIZE=${SHINYREPS_INSERTSIZE}" >> $output &&
            echo "SHINYREPS_GO_ENRICHMENT=${SHINYREPS_GO_ENRICHMENT}" >> $output &&
            echo "SHINYREPS_TRACKHUB_DONE=${SHINYREPS_TRACKHUB_DONE}" >> $output &&
            echo "SHINYREPS_TOOL_VERSIONS=${SHINYREPS_TOOL_VERSIONS}" >> $output &&
            echo "SHINYREPS_CUTADAPT_LOGS=${SHINYREPS_CUTADAPT_LOGS}" >> $output &&
            echo "SHINYREPS_GTF=${SHINYREPS_GTF}" >> $output &&
            echo "SHINYREPS_TARGET=${SHINYREPS_TARGET}" >> $output
        ""","shinyReports"
    }
}

