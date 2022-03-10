shinyReports = {
    doc title: "shinyReports",
        desc:  "creates the source code to compile the shiny and markdown reports",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Oliver Drechsel"

    output.dir = REPORTS

    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    produce("shinyReports.txt") {
        exec """
            ${PREAMBLE} &&

            if [ -e "${REPORTS}/variantreport.Rmd" ]; then
                echo 'variantreport.Rmd already exists. Older copy will be kept and not overwritten';
            else
                cp ${PIPELINE_ROOT}/tools/reports/shiny_dnaseq_reporting_tool/variantreport.Rmd ${REPORTS};
            fi &&

            cp ${PIPELINE_ROOT}/tools/reports/shiny_dnaseq_reporting_tool/variant.shinyrep.helpers.R ${REPORTS} &&
            cp ${PIPELINE_ROOT}/tools/reports/shiny_dnaseq_reporting_tool/styles.css ${REPORTS}                 &&

            PROJECT=\$(basename ${shinyReports_vars.project})                                 &&
            sed -i "2,2s/SHINYREPS_PROJECT/\${PROJECT}/" ${REPORTS}/variantreport.Rmd &&

            echo "SHINYREPS_PROJECT=${shinyReports_vars.project}" >  $output &&
            echo "SHINYREPS_PAIRED=${shinyReports_vars.paired}"   >> $output &&
            echo "SHINYREPS_LOG=$shinyReports_vars.log}"          >> $output &&
            echo "SHINYREPS_QC=${shinyReports_vars.qc}"           >> $output &&
            echo "SHINYREPS_RES=${shinyReports_vars.res}"         >> $output &&
            echo "SHINYREPS_PREFIX=${shinyReports_vars.prefix}"   >> $output &&
            echo "SHINYREPS_PLOTS_COLUMN=${shinyReports_vars.plots_column}"      >> $output &&
            echo "SHINYREPS_SORT_ALPHA=${shinyReports_vars.sort_alpha}"          >> $output &&
            echo "SHINYREPS_FASTQC_OUT=${shinyReports_vars.fastqc_out}"     >> $output &&
            echo "SHINYREPS_FASTQC_LOG=${shinyReports_vars.fastqc_log}"     >> $output &&
            echo "SHINYREPS_FASTQC_SUMMARIZED=${shinyReports_vars.fastqc_summarized}" >> $output &&
            echo "SHINYREPS_FASTQSCREEN_OUT=${shinyReports_vars.fastqscreen_out}" >> $output &&
            echo "SHINYREPS_FASTQSCREEN_LOG=${shinyReports_vars.fastqscreen_log}" >> $output &&
            echo "SHINYREPS_FASTQSCREEN_PERC=${shinyReports_vars.fastqscreen_perc}" >> $output &&
            echo "SHINYREPS_RUN_CUTADAPT=${shinyReports_vars.run_cutadapt}"  >> $output &&
 	    echo "SHINYREPS_CUTADAPT_STATS=${shinyReports_vars.cutadapt_stats}" >> $output &&
 	    echo "SHINYREPS_BQSR_STATS=${shinyReports_vars.bqsr_stats}" >> $output &&
            echo "SHINYREPS_BWA_LOG=${shinyReports_vars.bwa_log}"           >> $output &&
            echo "SHINYREPS_BWA_SUFFIX=${shinyReports_vars.bwa_suffix}"     >> $output &&
            echo "SHINYREPS_MARKDUPS_LOG=${shinyReports_vars.markdups_log}" >> $output &&
            echo "SHINYREPS_GATKhc_LOG=${shinyReports_vars.gatkhc_log}"           >> $output &&
            echo "SHINYREPS_GATKhc_SUFFIX=${shinyReports_vars.gatkhc_suffix}"     >> $output &&
            echo "SHINYREPS_GATKvariantcallingmetrics=${shinyReports_vars.variantcallingmetrics}" >> $output &&
            echo "SHINYREPS_GATKvarianteval=${shinyReports_vars.gatkvarianteval}" >> $output &&
            echo "SHINYREPS_GATKvarianteval_SUFFIX=${shinyReports_vars.gatkvarianteval_suffix}" >> $output &&
            echo "SHINYREPS_GATKhc_OUT=${shinyReports_vars.res_gatkhc}"     >> $output &&
            echo "SHINYREPS_RES_GATKhc_SUFFIX=${shinyReports_vars.res_gatkhc_suffix}" >> $output &&
            echo "SHINYREPS_SNPEFF_OUT=${shinyReports_vars.res_snpeff}"     >> $output &&
            echo "SHINYREPS_TOOL_VERSIONS=${shinyReports_vars.tool_versions}" >> $output &&
            echo "SHINYREPS_TARGET=${shinyReports_vars.target}" >> $output
        ""","shinyReports"
    }
}

