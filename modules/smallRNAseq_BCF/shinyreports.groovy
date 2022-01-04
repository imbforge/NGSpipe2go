shinyReports = {
    doc title: "shinyReports",
        desc:  "creates the source code to compile the shiny and markdown reports",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols, Anke Busch"

    output.dir = REPORTS

    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    produce("shinyReports.txt") {
        exec """
            ${PREAMBLE} &&

            cp ${PIPELINE_ROOT}/tools/reports/shiny_smallrnaseq_reporting_tool/smallRNA.shinyrep.helpers.R ${output.dir}   &&
            cp ${PIPELINE_ROOT}/tools/reports/shiny_smallrnaseq_reporting_tool/styles.css ${output.dir}              &&

            if [ -e "${output.dir}/smallRNAreport.Rmd" ]; then
                echo 'smallRNAreport.Rmd already exists. Older copy will be kept and not overwritten';
            else
                cp ${PIPELINE_ROOT}/tools/reports/shiny_smallrnaseq_reporting_tool/smallRNAreport.Rmd ${output.dir};
            fi &&

            if [ -e "${output.dir}/smallRNAreport.${shinyReports_vars.smallrna_type}.Rmd" ]; then
                echo 'smallRNAreport.${shinyReports_vars.smallrna_type}.Rmd already exists. Older copy will be kept and not overwritten';
            else
                cp ${PIPELINE_ROOT}/tools/reports/shiny_smallrnaseq_reporting_tool/smallRNAreport.type.Rmd ${output.dir}/smallRNAreport.${shinyReports_vars.smallrna_type}.Rmd;
            fi &&

            if [ ${shinyReports_vars.maturemirna} = "TRUE" ]; then
                if [ -e "${output.dir}/smallRNAreport.miRNAmature.Rmd" ]; then
                   echo 'smallRNAreport.miRNAmature.Rmd already exists. Older copy will be kept and not overwritten';
                else
                   cp ${PIPELINE_ROOT}/tools/reports/shiny_smallrnaseq_reporting_tool/smallRNAreport.miRNAmature.Rmd ${output.dir};
                fi;
            fi &&

            PROJECT=\$(basename ${shinyReports_vars.project})                            &&
            sed -i "2,2s/SHINYREPS_PROJECT/\${PROJECT}/" ${output.dir}/smallRNAreport.Rmd &&
            sed -i "2,2s/SHINYREPS_PROJECT/\${PROJECT}/" ${output.dir}/smallRNAreport.${shinyReports_vars.smallrna_type}.Rmd &&

            if [ ${shinyReports_vars.maturemirna} = "TRUE" ]; then
                sed -i "2,2s/SHINYREPS_PROJECT/\${PROJECT}/" ${output.dir}/smallRNAreport.miRNAmature.Rmd;
            fi &&
            
            echo "SHINYREPS_PROJECT=${shinyReports_vars.project}"                  > $output && 
            echo "SHINYREPS_ORG=${shinyReports_vars.org}"                         >> $output &&
            echo "SHINYREPS_PREFIX=${shinyReports_vars.prefix}"                   >> $output &&
            echo "SHINYREPS_PAIRED=${shinyReports_vars.paired}"                   >> $output &&
            echo "SHINYREPS_STRANDEDNESS=${shinyReports_vars.strandedness}"       >> $output &&
            echo "SHINYREPS_MINADAPTEROVERLAP=${shinyReports_vars.minadapteroverlap}" >> $output &&
            echo "SHINYREPS_MINREADLENGTH=${shinyReports_vars.minreadlength}"         >> $output &&
            echo "SHINYREPS_MAXREADLENGTH=${shinyReports_vars.maxreadlength}"         >> $output &&
            echo "SHINYREPS_CUTADAPT_STATS=${shinyReports_vars.cutadapt_stats}"       >> $output &&
            echo "SHINYREPS_QUALITYFILTER_LOG=${shinyReports_vars.qualityfilter_log}" >> $output &&
            echo "SHINYREPS_MINIMAL_QUAL=${shinyReports_vars.quality_min}"            >> $output &&
            echo "SHINYREPS_DEDUP_LOG=${shinyReports_vars.dedup_log}"                 >> $output &&
            echo "SHINYREPS_REMOVE_DUPLICATES=${shinyReports_vars.remove_dup}"        >> $output &&
            echo "SHINYREPS_BOWTIE_RES=${shinyReports_vars.bowtie_res}"               >> $output &&
            echo "SHINYREPS_BOWTIE_SUFFIX=${shinyReports_vars.bowtie_suffix}"         >> $output &&
            echo "SHINYREPS_FASTQC_OUT=${shinyReports_vars.fastqc_out}"               >> $output &&
            echo "SHINYREPS_FASTQC_SUMMARIZED=${shinyReports_vars.fastqc_summarized}" >> $output &&
            echo "SHINYREPS_FASTQSCREEN_OUT=${shinyReports_vars.fastqscreen_out}"     >> $output &&
            echo "SHINYREPS_FASTQSCREEN_PERC=${shinyReports_vars.fastqscreen_perc}"   >> $output &&
            echo "SHINYREPS_RNATYPES=${shinyReports_vars.rnatypes}"               >> $output &&
            echo "SHINYREPS_RNATYPES_SUFFIX=${shinyReports_vars.rnatypes_suffix}" >> $output &&
            echo "SHINYREPS_RNATYPES_CUTOFF=${shinyReports_vars.rnatypes_cutoff}" >> $output &&
            echo "SHINYREPS_SUBREAD=${shinyReports_vars.subread}"                 >> $output &&
            echo "SHINYREPS_SUBREAD_SUFFIX=${shinyReports_vars.subread_suffix}"   >> $output &&
            echo "SHINYREPS_SUBREAD_SUFFIX_TYPE=${shinyReports_vars.subread_suf_type}" >> $output &&
            echo "SHINYREPS_PLOTS_COLUMN=${shinyReports_vars.plots_column}"       >> $output &&
            echo "SHINYREPS_TOOL_VERSIONS=${shinyReports_vars.tool_versions}"     >> $output &&
            echo "SHINYREPS_TARGET=${shinyReports_vars.target}"                   >> $output &&
            echo "SHINYREPS_DE_DESEQ_DIR=${shinyReports_vars.de_deseq_dir}"       >> $output &&
            echo "SHINYREPS_DE_DESEQ_FILE=${shinyReports_vars.de_deseq_file}"     >> $output &&
            echo "SHINYREPS_DE_DESEQ_FDR=${shinyReports_vars.de_deseq_FDR}"       >> $output &&
            echo "SHINYREPS_DE_DESEQ_FC=${shinyReports_vars.de_deseq_FC}"         >> $output &&
            echo "SHINYREPS_GENETYPE=${shinyReports_vars.feature_type}"           >> $output &&
            echo "SHINYREPS_SMALLRNATYPE=${shinyReports_vars.smallrna_type}"      >> $output 
        ""","shinyReports"
    }
}

