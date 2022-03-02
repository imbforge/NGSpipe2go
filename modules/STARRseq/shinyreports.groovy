shinyReports = {
    doc title: "shinyReports",
        desc:  "creates the source code to compile the shiny and markdown reports",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols"

    output.dir = REPORTS

    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    produce("shinyReports.txt") {
        exec """
            ${PREAMBLE} &&

            cp ${PIPELINE_ROOT}/tools/reports/shiny_starrseq_reporting_tool/styles.css ${REPORTS}              &&

            if [ ${PIPELINE} == "CapSTARRseq" ]; then
                cp ${PIPELINE_ROOT}/tools/reports/shiny_starrseq_reporting_tool/CapSTARR.shinyrep.helpers.R ${REPORTS};
                if [ -e "${REPORTS}/CapSTARRreport.Rmd" ]; then
                    echo 'CapSTARRreport.Rmd already exists. Older copy will be kept and not overwritten';
                else
                    cp ${PIPELINE_ROOT}/tools/reports/shiny_starrseq_reporting_tool/CapSTARRreport.Rmd ${REPORTS};
                fi
            else
                cp ${PIPELINE_ROOT}/tools/reports/shiny_starrseq_reporting_tool/STARR.shinyrep.helpers.R ${REPORTS};
                if [ -e "${REPORTS}/STARRreport.Rmd" ]; then
                    echo 'STARRreport.Rmd already exists. Older copy will be kept and not overwritten';
                else
                    cp ${PIPELINE_ROOT}/tools/reports/shiny_starrseq_reporting_tool/STARRreport.Rmd ${REPORTS};
                fi
            fi &&


            PROJECT=\$(basename ${shinyReports_vars.project})                    &&
            if [ ${PIPELINE} == "CapSTARRseq" ]; then
                sed -i "2,2s/SHINYREPS_PROJECT/\${PROJECT}/" ${REPORTS}/CapSTARRreport.Rmd;
            else
                sed -i "2,2s/SHINYREPS_PROJECT/\${PROJECT}/" ${REPORTS}/STARRreport.Rmd;
            fi &&

            echo "SHINYREPS_PROJECT=${shinyReports_vars.project}" >  $output &&
            echo "SHINYREPS_ORG=${shinyReports_vars.org}"         >> $output &&
            echo "SHINYREPS_DB=${shinyReports_vars.db}"           >> $output &&
            echo "SHINYREPS_LOG=${shinyReports_vars.log}"         >> $output &&
            echo "SHINYREPS_PAIRED=${shinyReports_vars.paired}"   >> $output &&
            echo "SHINYREPS_QC=${shinyReports_vars.qc}"           >> $output &&
            echo "SHINYREPS_RES=${shinyReports_vars.res}"         >> $output &&
            echo "SHINYREPS_PREFIX=${shinyReports_vars.prefix}"   >> $output &&
            echo "SHINYREPS_STRANDEDNESS=${shinyReports_vars.strandedness}"  >> $output &&
            echo "SHINYREPS_RUN_CUTADAPT=${shinyReports_vars.run_cutadapt}"  >> $output &&
            echo "SHINYREPS_RUN_PEAK_ANNOTATION=${shinyReports_vars.run_peakanno}"  >> $output &&
            echo "SHINYREPS_RUN_DIFFBIND=${shinyReports_vars.run_diffbind}"  >> $output &&
            echo "SHINYREPS_RUN_ENRICHMENT=${shinyReports_vars.run_enrich}"  >> $output &&
            echo "SHINYREPS_CAPSTARRSEQ_DIFFEXP=${shinyReports_vars.capstarrseq_diffexp}"  >> $output &&
            echo "SHINYREPS_GTF=${shinyReports_vars.gtf}"         >> $output &&
            echo "SHINYREPS_FASTQC_OUT=${shinyReports_vars.fastqc_out}"     >> $output &&
            echo "SHINYREPS_FASTQC_LOG=${shinyReports_vars.fastqc_log}"     >> $output &&
            echo "SHINYREPS_FASTQC_SUMMARIZED=${shinyReports_vars.fastqc_summarized}" >> $output &&
            echo "SHINYREPS_FASTQSCREEN_OUT=${shinyReports_vars.fastqscreen_out}" >> $output &&
            echo "SHINYREPS_FASTQSCREEN_LOG=${shinyReports_vars.fastqscreen_log}" >> $output &&
            echo "SHINYREPS_FASTQSCREEN_PERC=${shinyReports_vars.fastqscreen_perc}" >> $output &&
            echo "SHINYREPS_CUTADAPT_STATS=${shinyReports_vars.cutadapt_stats}" >> $output &&
            echo "SHINYREPS_BOWTIE_LOG=${shinyReports_vars.bowtie_log}"      >> $output &&
            echo "SHINYREPS_BAMINDEX_LOG=${shinyReports_vars.bamindex_log}"     >> $output &&
            echo "SHINYREPS_INSERTSIZE=${shinyReports_vars.insertsize}" >> $output &&
            echo "SHINYREPS_EXTEND_LOG=${shinyReports_vars.extend_log}"      >> $output &&
            echo "SHINYREPS_DUPRADAR_LOG=${shinyReports_vars.dupradar_log}" >> $output &&
            echo "SHINYREPS_DE_DESEQ=${shinyReports_vars.de_deseq}"       >> $output &&
            echo "SHINYREPS_DE_DESEQ_FDR=${shinyReports_vars.de_deseq_FDR}"       >> $output &&
            echo "SHINYREPS_DE_DESEQ_FC=${shinyReports_vars.de_deseq_FC}"       >> $output &&
            echo "SHINYREPS_SUBREAD=${shinyReports_vars.subread}"         >> $output &&
            echo "SHINYREPS_SUBREAD_SUFFIX=${shinyReports_vars.subread_suffix}"  >> $output &&
            echo "SHINYREPS_SUBREAD_LOG=${shinyReports_vars.subread_log}"        >> $output &&
            echo "SHINYREPS_BAM2BW_LOG=${shinyReports_vars.bam2bw_log}"          >> $output &&
            echo "SHINYREPS_MARKDUPS_LOG=${shinyReports_vars.markdups_log}"      >> $output &&
            echo "SHINYREPS_DESEQ_LOGS=${shinyReports_vars.deseq_logs}"          >> $output &&
            echo "SHINYREPS_PLOTS_COLUMN=${shinyReports_vars.plots_column}"      >> $output &&
            echo "SHINYREPS_SORT_ALPHA=${shinyReports_vars.sort_alpha}"          >> $output &&
            echo "SHINYREPS_IPSTRENGTH=${shinyReports_vars.ipstrength}"      >> $output &&
            echo "SHINYREPS_IPSTRENGTH_LOG=${shinyReports_vars.ipstrength_log}"      >> $output &&
            echo "SHINYREPS_PBC=${shinyReports_vars.pbc}"         >> $output &&
            echo "SHINYREPS_PHANTOMPEAK=${shinyReports_vars.phantompeak}"    >> $output &&
            echo "SHINYREPS_PHANTOM_LOG=${shinyReports_vars.phantom_log}"    >> $output &&
            echo "SHINYREPS_MACS2=${shinyReports_vars.macs2}"     >> $output &&
            echo "SHINYREPS_MACS2_LOG=${shinyReports_vars.macs2_log}"         >> $output &&
            echo "SHINYREPS_BLACKLIST_FILTER=${shinyReports_vars.blacklist_filter}" >> $output &&
            echo "SHINYREPS_UPSETPLOT=${shinyReports_vars.upset}" >> $output &&
            echo "SHINYREPS_PEAK_ANNOTATION=${shinyReports_vars.peak_annotation}" >> $output &&
            echo "SHINYREPS_GREAT=${shinyReports_vars.great}" >> $output &&
            echo "SHINYREPS_DIFFBIND=${shinyReports_vars.diffbind}" >> $output           &&
            echo "SHINYREPS_TRACKHUB_DONE=${shinyReports_vars.trackhub_done}" >> $output &&
            echo "SHINYREPS_TOOL_VERSIONS=${shinyReports_vars.tool_versions}" >> $output &&
            echo "SHINYREPS_TARGET=${shinyReports_vars.target}" >> $output
        ""","shinyReports"
    }
}

