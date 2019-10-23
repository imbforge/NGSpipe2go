// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/shinyreports.vars.groovy"

shinyReports = {
    doc title: "shinyReports",
        desc:  "creates the source code to compile the shiny and markdown reports",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols, Anke Busch"

    output.dir = REPORTS

    def PREAMBLE = get_preamble("shinyReports")

    produce("shinyReports.txt") {
        exec """
            ${PREAMBLE} &&

            cp ${PIPELINE_ROOT}/tools/reports/shiny_smallrnaseq_reporting_tool/smallRNA.shinyrep.helpers.R ${REPORTS}   &&
            cp ${PIPELINE_ROOT}/tools/reports/shiny_smallrnaseq_reporting_tool/styles.css ${REPORTS}              &&

            if [ -e "${REPORTS}/smallRNAreport.Rmd" ]; then
                echo 'smallRNAreport.Rmd already exists. Older copy will be kept and not overwritten';
            else
                cp ${PIPELINE_ROOT}/tools/reports/shiny_smallrnaseq_reporting_tool/smallRNAreport.Rmd ${REPORTS};
            fi &&

            PROJECT=\$(basename ${shinyReports_vars.project})                            &&
            sed -i "2,2s/SHINYREPS_PROJECT/\${PROJECT}/" ${REPORTS}/smallRNAreport.Rmd &&

            echo "SHINYREPS_PROJECT=${shinyReports_vars.project}"                  > $output &&
            echo "SHINYREPS_MINADAPTEROVERLAP=${shinyReports_vars.minadapteroverlap}" >> $output &&
            echo "SHINYREPS_MINREADLENGTH=${shinyReports_vars.minreadlength}"     >> $output &&
            echo "SHINYREPS_PREFIX=${shinyReports_vars.prefix}"                   >> $output &&
            echo "SHINYREPS_BOWTIE_RES=${shinyReports_vars.bowtie_res}"           >> $output &&
            echo "SHINYREPS_BOWTIE_LOG=${shinyReports_vars.bowtie_log}"           >> $output &&
            echo "SHINYREPS_BOWTIE_SUFFIX=${shinyReports_vars.bowtie_suffix}"     >> $output &&
            echo "SHINYREPS_BOWTIE_PLOT_DIR=${shinyReports_vars.bowtie_plot_dir}" >> $output &&
            echo "SHINYREPS_STAR_LOG=${shinyReports_vars.star_log}"               >> $output &&
            echo "SHINYREPS_STAR_SUFFIX=${shinyReports_vars.star_suffix}"         >> $output &&
            echo "SHINYREPS_STARparms_SUFFIX=${shinyReports_vars.starparms_suffix}"   >> $output &&
            echo "SHINYREPS_STAR_PLOT_DIR=${shinyReports_vars.star_plot_dir}"     >> $output &&
            echo "SHINYREPS_CUTADAPT_PLOT_DIR=${shinyReports_vars.cutadapt_plot_dir}" >> $output &&
            echo "SHINYREPS_CUTADAPT_LOG=${shinyReports_vars.cutadapt_log}"       >> $output &&
            echo "SHINYREPS_QUALITYFILTER_PLOT_DIR=${shinyReports_vars.qualityfilter_plot_dir}" >> $output &&
            echo "SHINYREPS_QUALITYFILTER_LOG=${shinyReports_vars.qualityfilter_log}" >> $output &&
            echo "SHINYREPS_DEDUP_PLOT_DIR=${shinyReports_vars.dedup_plot_dir}"   >> $output &&
            echo "SHINYREPS_DEDUP_LOG=${shinyReports_vars.dedup_log}"             >> $output &&
            echo "SHINYREPS_RAWFILTERSUMMARY_PLOT_DIR=${shinyReports_vars.rawfiltersummary_plot_dir}" >> $output &&
            echo "SHINYREPS_FASTQC_OUT=${shinyReports_vars.fastqc_out}"           >> $output &&
            echo "SHINYREPS_FASTQC_LOG=${shinyReports_vars.fastqc_log}"           >> $output &&
            echo "SHINYREPS_FASTQSCREEN_OUT=${shinyReports_vars.fastqscreen_out}" >> $output &&
            echo "SHINYREPS_FASTQSCREEN_LOG=${shinyReports_vars.fastqscreen_log}" >> $output &&
            echo "SHINYREPS_BAMINDEX_LOG=${shinyReports_vars.bamindex_log}"       >> $output &&
            echo "SHINYREPS_RNATYPES_LOG=${shinyReports_vars.rnatypes_log}"       >> $output &&
            echo "SHINYREPS_RNATYPES=${shinyReports_vars.rnatypes}"               >> $output &&
            echo "SHINYREPS_RNATYPES_SUFFIX=${shinyReports_vars.rnatypes_suffix}" >> $output &&
            echo "SHINYREPS_RNATYPES_CUTOFF=${shinyReports_vars.rnatypes_cutoff}" >> $output &&
            echo "SHINYREPS_SUBREAD=${shinyReports_vars.subread}"                 >> $output &&
            echo "SHINYREPS_SUBREAD_SUFFIX=${shinyReports_vars.subread_suffix}"   >> $output &&
            echo "SHINYREPS_SUBREAD_LOG=${shinyReports_vars.subread_log}"         >> $output &&
            echo "SHINYREPS_PLOTS_COLUMN=${shinyReports_vars.plots_column}"       >> $output &&
            echo "SHINYREPS_TOOL_VERSIONS=${shinyReports_vars.tool_versions}"     >> $output 
        ""","shinyReports"
    }
}

