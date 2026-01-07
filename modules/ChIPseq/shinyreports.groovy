shinyReports = {
    doc title: "shinyReports",
        desc:  "creates the source code to compile the shiny and markdown reports",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Sergi Sayols"

    output.dir = REPORTS

    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    produce("shinyReports.txt") {
        exec """
            ${PREAMBLE} &&

            cp ${PIPELINE_ROOT}/tools/reports/shiny_chipseq_reporting_tool/ChIP.shinyrep.helpers.R ${REPORTS} &&
            cp ${PIPELINE_ROOT}/tools/reports/shiny_chipseq_reporting_tool/styles.css ${REPORTS}              &&

            if [ -e "${REPORTS}/ChIPreport.Rmd" ]; then
                echo 'ChIPreport.Rmd already exists. Older copy will be kept and not overwritten';
            else
                cp ${PIPELINE_ROOT}/tools/reports/shiny_chipseq_reporting_tool/ChIPreport.Rmd ${REPORTS};
            fi &&

            PROJECT=\$(basename ${shinyReports_vars.project})                              &&
            sed -i "2,2s/SHINYREPS_PROJECT/\${PROJECT}/" ${REPORTS}/ChIPreport.Rmd &&

            echo "SHINYREPS_PROJECT=${shinyReports_vars.project}" >  $output &&
            echo "SHINYREPS_LOG=${shinyReports_vars.log}"         >> $output &&
            echo "SHINYREPS_QC=${shinyReports_vars.qc}"           >> $output &&
            echo "SHINYREPS_RES=${shinyReports_vars.res}"         >> $output &&
            echo "SHINYREPS_TARGET=${shinyReports_vars.target}"   >> $output &&
            echo "SHINYREPS_PAIRED=${shinyReports_vars.paired}"      >> $output &&
            echo "SHINYREPS_RUN_CUTADAPT=${shinyReports_vars.run_cutadapt}"  >> $output &&
            echo "SHINYREPS_RUN_PEAK_ANNOTATION=${shinyReports_vars.run_peakanno}"  >> $output &&
            echo "SHINYREPS_RUN_DIFFBIND=${shinyReports_vars.run_diffbind}"  >> $output &&
            echo "SHINYREPS_RUN_ENRICHMENT=${shinyReports_vars.run_enrich}"  >> $output &&
 	        echo "SHINYREPS_CUTADAPT_STATS=${shinyReports_vars.cutadapt_stats}" >> $output &&
 	        echo "SHINYREPS_SAMTOOLSCOV_OUT=${shinyReports_vars.samtoolscov_out}" >> $output &&
 	        echo "SHINYREPS_BOWTIE_LOG=${shinyReports_vars.bowtie_log}"      >> $output &&
            echo "SHINYREPS_BAMINDEX_LOG=${shinyReports_vars.bamindex_log}"  >> $output &&
            echo "SHINYREPS_MARKDUPS_LOG=${shinyReports_vars.markdups_log}"  >> $output &&
            echo "SHINYREPS_EXTEND_LOG=${shinyReports_vars.extend_log}"      >> $output &&
            echo "SHINYREPS_FASTQC=${shinyReports_vars.fastqc}"   >> $output &&
            echo "SHINYREPS_FASTQC_LOG=${shinyReports_vars.fastqc_log}"   >> $output &&
            echo "SHINYREPS_FASTQC_SUMMARIZED=${shinyReports_vars.fastqc_summarized}" >> $output &&
            echo "SHINYREPS_FASTQSCREEN_OUT=${shinyReports_vars.fastqscreen_out}" >> $output &&
            echo "SHINYREPS_FASTQSCREEN_LOG=${shinyReports_vars.fastqscreen_log}" >> $output &&
            echo "SHINYREPS_FASTQSCREEN_PERC=${shinyReports_vars.fastqscreen_perc}" >> $output &&
            echo "SHINYREPS_INSERTSIZE=${shinyReports_vars.insertsize}" >> $output &&
            echo "SHINYREPS_IPSTRENGTH=${shinyReports_vars.ipstrength}"      >> $output &&
            echo "SHINYREPS_IPSTRENGTH_LOG=${shinyReports_vars.ipstrength_log}"      >> $output &&
            echo "SHINYREPS_PBC=${shinyReports_vars.pbc}"         >> $output &&
            echo "SHINYREPS_PHANTOMPEAK=${shinyReports_vars.phantompeak}"    >> $output &&
            echo "SHINYREPS_PHANTOM_LOG=${shinyReports_vars.phantom_log}"    >> $output &&
            echo "SHINYREPS_MACS2=${shinyReports_vars.macs2}"     >> $output &&
            echo "SHINYREPS_MACS2_LOG=${shinyReports_vars.macs2_log}"         >> $output &&
            echo "SHINYREPS_EXCLUDEDREGIONS_FILTER=${shinyReports_vars.excludedRegions_filter}" >> $output &&
            echo "SHINYREPS_UPSETPLOT=${shinyReports_vars.upset}" >> $output &&
            echo "SHINYREPS_PREFIX=${shinyReports_vars.prefix}"   >> $output &&
            echo "SHINYREPS_PLOTS_COLUMN=${shinyReports_vars.plots_column}" >> $output &&
            echo "SHINYREPS_PEAK_ANNOTATION=${shinyReports_vars.peak_annotation}" >> $output &&
            echo "SHINYREPS_DB=${shinyReports_vars.db}"   >> $output &&
            echo "SHINYREPS_GREAT=${shinyReports_vars.great}" >> $output &&
            echo "SHINYREPS_DIFFBIND=${shinyReports_vars.diffbind}" >> $output           &&
            echo "SHINYREPS_TRACKHUB_DONE=${shinyReports_vars.trackhub_done}" >> $output &&
            echo "SHINYREPS_TOOL_VERSIONS=${shinyReports_vars.tool_versions}" >> $output
        ""","shinyReports"
    }
}

