//rule for task shinyReports from catalog smallRNAseq_BCF, version 0.1
//desc: creates the source code to compile the shiny and markdown reports
shinyReports = {
	doc title: "shinyReports",
		desc:  "creates the source code to compile the shiny and markdown reports",
		constraints: "",
		bpipe_version: "tested with bpipe 0.9.8.7",
		author: "Sergi Sayols, Anke Busch"
	
	output.dir = REPORTS
	
	produce("shinyReports.txt") {
		exec """
			
			cp ${MODULE_FOLDER}/../tools/reports/shiny_smallrnaseq_reporting_tool/smallRNA.shinyrep.helpers.R ${REPORTS}   &&
			cp ${MODULE_FOLDER}/../tools/reports/shiny_smallrnaseq_reporting_tool/styles.css ${REPORTS}              &&
			
			if [ -e "${REPORTS}/smallRNAreport.Rmd" ]; then
				echo 'smallRNAreport.Rmd already exists. Older copy will be kept and not overwritten';
			else
				cp ${MODULE_FOLDER}/../tools/reports/shiny_smallrnaseq_reporting_tool/smallRNAreport.Rmd ${REPORTS};
			fi &&
			
			PROJECT=\$(basename ${SHINYREPS_PROJECT})                            &&
			sed -i "2,2s/SHINYREPS_PROJECT/\${PROJECT}/" ${REPORTS}/smallRNAreport.Rmd &&
			
			echo "SHINYREPS_PROJECT=${SHINYREPS_PROJECT}" >  $output &&
			echo "SHINYREPS_MINADAPTEROVERLAP=${SHINYREPS_MINADAPTEROVERLAP}" >> $output &&
			echo "SHINYREPS_MINREADLENGTH=${SHINYREPS_MINREADLENGTH}" >> $output &&
			echo "SHINYREPS_PREFIX=${SHINYREPS_PREFIX}"   >> $output &&
			echo "SHINYREPS_BOWTIE_RES=${SHINYREPS_BOWTIE_RES}" >> $output &&
			echo "SHINYREPS_BOWTIE_LOG=${SHINYREPS_BOWTIE_LOG}" >> $output &&
			echo "SHINYREPS_BOWTIE_SUFFIX=${SHINYREPS_BOWTIE_SUFFIX}" >> $output &&
			echo "SHINYREPS_BOWTIE_PLOT_DIR=${SHINYREPS_BOWTIE_PLOT_DIR}" >> $output &&
			echo "SHINYREPS_STAR_LOG=${SHINYREPS_STAR_LOG}"       >> $output &&
                        echo "SHINYREPS_STAR_SUFFIX=${SHINYREPS_STAR_SUFFIX}" >> $output &&
                        echo "SHINYREPS_STARparms_SUFFIX=${SHINYREPS_STARparms_SUFFIX}" >> $output &&
			echo "SHINYREPS_STAR_PLOT_DIR=${SHINYREPS_STAR_PLOT_DIR}" >> $output &&
			echo "SHINYREPS_CUTADAPT_PLOT_DIR=${SHINYREPS_CUTADAPT_PLOT_DIR}" >> $output &&
			echo "SHINYREPS_CUTADAPT_LOG=${SHINYREPS_CUTADAPT_LOG}" >> $output &&
			echo "SHINYREPS_QUALITYFILTER_PLOT_DIR=${SHINYREPS_QUALITYFILTER_PLOT_DIR}" >> $output &&
			echo "SHINYREPS_QUALITYFILTER_LOG=${SHINYREPS_QUALITYFILTER_LOG}" >> $output &&
			echo "SHINYREPS_DEDUP_PLOT_DIR=${SHINYREPS_DEDUP_PLOT_DIR}" >> $output &&
			echo "SHINYREPS_DEDUP_LOG=${SHINYREPS_DEDUP_LOG}" >> $output &&
			echo "SHINYREPS_RAWFILTERSUMMARY_PLOT_DIR=${SHINYREPS_RAWFILTERSUMMARY_PLOT_DIR}" >> $output &&
			echo "SHINYREPS_FASTQC_OUT=${SHINYREPS_FASTQC_OUT}"     >> $output &&
			echo "SHINYREPS_FASTQC_LOG=${SHINYREPS_FASTQC_LOG}"     >> $output &&
			echo "SHINYREPS_FASTQSCREEN_OUT=${SHINYREPS_FASTQSCREEN_OUT}"     >> $output &&
			echo "SHINYREPS_FASTQSCREEN_LOG=${SHINYREPS_FASTQSCREEN_LOG}"     >> $output &&
			echo "SHINYREPS_BAMINDEX_LOG=${SHINYREPS_BAMINDEX_LOG}"     >> $output &&
			echo "SHINYREPS_RNATYPES_LOG=${SHINYREPS_RNATYPES_LOG}" >> $output &&
			echo "SHINYREPS_RNATYPES=${SHINYREPS_RNATYPES}" >> $output &&
			echo "SHINYREPS_RNATYPES_SUFFIX=${SHINYREPS_RNATYPES_SUFFIX}" >> $output &&
			echo "SHINYREPS_RNATYPES_CUTOFF=${SHINYREPS_RNATYPES_CUTOFF}" >> $output &&
			echo "SHINYREPS_SUBREAD=${SHINYREPS_SUBREAD}"         >> $output &&
			echo "SHINYREPS_SUBREAD_SUFFIX=${SHINYREPS_SUBREAD_SUFFIX}"  >> $output &&
			echo "SHINYREPS_SUBREAD_LOG=${SHINYREPS_SUBREAD_LOG}"        >> $output &&
			echo "SHINYREPS_PLOTS_COLUMN=${SHINYREPS_PLOTS_COLUMN}"      >> $output 
		""","shinyReports"
	}
}

