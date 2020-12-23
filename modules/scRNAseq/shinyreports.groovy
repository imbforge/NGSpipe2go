shinyReports = {
    doc title: "shinyReports",
        desc:  "creates the source code to compile the shiny and markdown reports",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols, modified by Frank Rühle"

    output.dir = REPORTS

    def PREAMBLE = get_preamble("shinyReports")

    produce("shinyReports.txt") {
        exec """
            ${PREAMBLE} &&

            cp ${PIPELINE_ROOT}/tools/reports/shiny_scrnaseq_reporting_tool/server.R ${REPORTS}                &&
            cp ${PIPELINE_ROOT}/tools/reports/shiny_scrnaseq_reporting_tool/ui.R ${REPORTS}                    &&
            cp ${PIPELINE_ROOT}/tools/reports/shiny_scrnaseq_reporting_tool/sc.shinyrep.helpers.R ${REPORTS}   &&
            cp ${PIPELINE_ROOT}/tools/reports/shiny_scrnaseq_reporting_tool/bustard.pl ${REPORTS}              &&
            cp ${PIPELINE_ROOT}/tools/reports/shiny_scrnaseq_reporting_tool/BustardSummary.toMD.xsl ${REPORTS} &&
            cp ${PIPELINE_ROOT}/tools/reports/shiny_scrnaseq_reporting_tool/styles.css ${REPORTS}              &&
            cp ${PIPELINE_ROOT}/tools/reports/shiny_scrnaseq_reporting_tool/app.R ${REPORTS}                   &&

            if [ -e "${REPORTS}/sc.report.Rmd" ]; then
                echo 'sc.report.Rmd already exists. Older copy will be kept and not overwritten';
            else
                cp ${PIPELINE_ROOT}/tools/reports/shiny_scrnaseq_reporting_tool/sc.report.Rmd ${REPORTS};
            fi &&

            PROJECT=\$(basename ${shinyReports_vars.project})                     &&
            sed -i "2,2s/SHINYREPS_PROJECT/\${PROJECT}/" ${REPORTS}/sc.report.Rmd &&

            echo "SHINYREPS_PROJECT=${shinyReports_vars.project}" >  $output &&
            echo "SHINYREPS_ORG=${shinyReports_vars.org}"         >> $output &&
            echo "SHINYREPS_DB=${shinyReports_vars.db}"           >> $output &&
            echo "SHINYREPS_LOG=${shinyReports_vars.log}"         >> $output &&
            echo "SHINYREPS_PAIRED=${shinyReports_vars.paired}"   >> $output &&
            echo "SHINYREPS_QC=${shinyReports_vars.qc}"           >> $output &&
            echo "SHINYREPS_RES=${shinyReports_vars.res}"         >> $output &&
            echo "SHINYREPS_CELLRANGERAGGR_ID=${shinyReports_vars.cellranger_aggr_id}"  >> $output &&
            echo "SHINYREPS_PREFIX=${shinyReports_vars.prefix}"   >> $output &&
            echo "SHINYREPS_STRANDEDNESS=${shinyReports_vars.strandedness}" >> $output &&
            echo "SHINYREPS_STAR_LOG=${shinyReports_vars.star_log}"       >> $output &&
            echo "SHINYREPS_STAR_SUFFIX=${shinyReports_vars.star_suffix}" >> $output &&
            echo "SHINYREPS_STARparms_SUFFIX=${shinyReports_vars.starparms_suffix}" >> $output &&
            echo "SHINYREPS_FASTQC_OUT=${shinyReports_vars.fastqc_out}"     >> $output &&
            echo "SHINYREPS_FASTQC_LOG=${shinyReports_vars.fastqc_log}"     >> $output &&
            echo "SHINYREPS_FASTQC_SUMMARIZED=${shinyReports_vars.fastqc_summarized}" >> $output &&
            echo "SHINYREPS_FASTQSCREEN_OUT=${shinyReports_vars.fastqscreen_out}" >> $output &&
            echo "SHINYREPS_FASTQSCREEN_LOG=${shinyReports_vars.fastqscreen_log}" >> $output &&
            echo "SHINYREPS_BAMINDEX_LOG=${shinyReports_vars.bamindex_log}" >> $output &&
            echo "SHINYREPS_DUPRADAR_LOG=${shinyReports_vars.dupradar_log}" >> $output &&
            echo "SHINYREPS_RNATYPES_LOG=${shinyReports_vars.rnatypes_log}" >> $output &&
            echo "SHINYREPS_RNATYPES=${shinyReports_vars.rnatypes}" >> $output &&
            echo "SHINYREPS_RNATYPES_SUFFIX=${shinyReports_vars.rnatypes_suffix}" >> $output &&
            echo "SHINYREPS_GENEBODYCOV_LOG=${shinyReports_vars.genebodycov_log}" >> $output &&
            echo "SHINYREPS_BUSTARD=${shinyReports_vars.bustard}"         >> $output &&
            echo "SHINYREPS_SUBREAD=${shinyReports_vars.subread}"         >> $output &&
            echo "SHINYREPS_SUBREAD_SUFFIX=${shinyReports_vars.subread_suffix}" >> $output &&
            echo "SHINYREPS_SUBREAD_LOG=${shinyReports_vars.subread_log}"       >> $output &&
            echo "SHINYREPS_UMICOUNT=${shinyReports_vars.umicount}"             >> $output &&
            echo "SHINYREPS_UMICOUNT_LOG=${shinyReports_vars.umicount_log}"     >> $output &&
            echo "SHINYREPS_BAM2BW_LOG=${shinyReports_vars.bam2bw_log}"         >> $output &&
            echo "SHINYREPS_MARKDUPS_LOG=${shinyReports_vars.markdups_log}"     >> $output &&
            echo "SHINYREPS_PLOTS_COLUMN=${shinyReports_vars.plots_column}"     >> $output &&
            echo "SHINYREPS_INFEREXPERIMENT_LOGS=${shinyReports_vars.inferexperiment_logs}" >> $output &&
            echo "SHINYREPS_QUALIMAP_LOGS=${shinyReports_vars.qualimap_logs}" >> $output &&
            echo "SHINYREPS_INSERTSIZE=${shinyReports_vars.insertsize}" >> $output &&
            echo "SHINYREPS_TRACKHUB_DONE=${shinyReports_vars.trackhub_done}" >> $output &&
            echo "SHINYREPS_TOOL_VERSIONS=${shinyReports_vars.tool_versions}" >> $output &&
            echo "SHINYREPS_RUN_CUTADAPT=${shinyReports_vars.run_cutadapt}"  >> $output &&
            echo "SHINYREPS_CUTADAPT_LOGS=${shinyReports_vars.cutadapt_logs}" >> $output &&
            echo "SHINYREPS_GTF=${shinyReports_vars.gtf}"         >> $output &&
            echo "SHINYREPS_TARGET=${shinyReports_vars.target}"   >> $output &&
            echo "SHINYREPS_MTGENES=${shinyReports_vars.mtgenes}" >> $output &&
            echo "SHINYREPS_SAMPLEPATTERN1=${shinyReports_vars.samplepattern1}" >> $output &&
            echo "SHINYREPS_SAMPLEPATTERN2=${shinyReports_vars.samplepattern2}" >> $output &&
            echo "SHINYREPS_COLORBYFACTOR=${shinyReports_vars.colorByFactor}" >> $output &&
            echo "SHINYREPS_MAXNO=${shinyReports_vars.maxno}" >> $output &&
            echo "SHINYREPS_SEQTYPE=${shinyReports_vars.seqtype}" >> $output &&
            echo "SHINYREPS_TYPE_OF_THRESHOLD=${shinyReports_vars.type_of_threshold}" >> $output &&
            echo "SHINYREPS_THRESHOLD_TOTAL_COUNTS_MIN=${shinyReports_vars.threshold_total_counts_min}" >> $output &&
            echo "SHINYREPS_THRESHOLD_TOTAL_COUNTS_MAX=${shinyReports_vars.threshold_total_counts_max}" >> $output &&
            echo "SHINYREPS_THRESHOLD_TOTAL_FEATURES_BY_COUNTS_ENDOGENOUS=${shinyReports_vars.threshold_total_features_by_counts_endogenous}" >> $output &&
            echo "SHINYREPS_THRESHOLD_PCT_COUNTS_MT=${shinyReports_vars.threshold_pct_counts_Mt}" >> $output &&
            echo "SHINYREPS_NMADS=${shinyReports_vars.nmads}" >> $output
        ""","shinyReports"
    }
}

