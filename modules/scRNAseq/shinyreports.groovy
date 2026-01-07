shinyReports = {
    doc title: "shinyReports",
        desc:  "creates the source code to compile the shiny and markdown reports",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols, modified by Frank Rühle"

    output.dir = REPORTS
    println "${REPORTS}/${shinyReports_vars.report}"
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    produce("shinyReports.txt") {
        exec """
            ${PREAMBLE} &&

            cp ${PIPELINE_ROOT}/tools/reports/shiny_scrnaseq_reporting_tool/server.R ${REPORTS}                &&
            cp ${PIPELINE_ROOT}/tools/reports/shiny_scrnaseq_reporting_tool/ui.R ${REPORTS}                    &&
            cp ${PIPELINE_ROOT}/tools/reports/shiny_scrnaseq_reporting_tool/sc.shinyrep.helpers.R ${REPORTS}   &&
            cp ${PIPELINE_ROOT}/tools/reports/shiny_scrnaseq_reporting_tool/styles.css ${REPORTS}              &&
            cp ${PIPELINE_ROOT}/tools/reports/shiny_scrnaseq_reporting_tool/app.R ${REPORTS}                   &&

            if [ -e "${REPORTS}/${shinyReports_vars.report}" ]; then
                echo "${shinyReports_vars.report} already exists. Older copy will be kept and not overwritten";
            else
                cp ${PIPELINE_ROOT}/tools/reports/shiny_scrnaseq_reporting_tool/${shinyReports_vars.report} ${REPORTS};
                if [ "${shinyReports_vars.seqtype}" == "tenXmultiome" ]; then
                    cp ${PIPELINE_ROOT}/tools/reports/shiny_scrnaseq_reporting_tool/sc.report.Rmd ${REPORTS};
                    cp ${PIPELINE_ROOT}/tools/reports/shiny_scrnaseq_reporting_tool/scatac.report.Rmd ${REPORTS};
                fi
            fi &&

            PROJECT=\$(basename ${shinyReports_vars.project})                     &&
            sed -i "2,2s/SHINYREPS_PROJECT/\${PROJECT}/" ${REPORTS}/${shinyReports_vars.report} &&

            echo "SHINYREPS_PROJECT=${shinyReports_vars.project}" >  $output &&
            echo "SHINYREPS_SEQTYPE=${shinyReports_vars.seqtype}" >> $output &&
            echo "SHINYREPS_GTF=${shinyReports_vars.gtf}"         >> $output &&
            echo "SHINYREPS_TARGET=${shinyReports_vars.target}"   >> $output &&
            echo "SHINYREPS_CONTRASTS=${shinyReports_vars.contrasts}"   >> $output &&
            echo "SHINYREPS_ORG=${shinyReports_vars.org}"         >> $output &&
            echo "SHINYREPS_DB=${shinyReports_vars.db}"           >> $output &&
            echo "SHINYREPS_LOG=${shinyReports_vars.log}"         >> $output &&
            echo "SHINYREPS_PAIRED=${shinyReports_vars.paired}"   >> $output &&
            echo "SHINYREPS_QC=${shinyReports_vars.qc}"           >> $output &&
            echo "SHINYREPS_RES=${shinyReports_vars.res}"         >> $output &&
            echo "SHINYREPS_PREFIX=${shinyReports_vars.prefix}"   >> $output &&
            echo "SHINYREPS_STRANDEDNESS=${shinyReports_vars.strandedness}" >> $output &&
            echo "SHINYREPS_STAR_LOG=${shinyReports_vars.star_log}"       >> $output &&
            echo "SHINYREPS_STAR_SUFFIX=${shinyReports_vars.star_suffix}" >> $output &&
            echo "SHINYREPS_STARparms_SUFFIX=${shinyReports_vars.starparms_suffix}" >> $output &&
            echo "SHINYREPS_FASTQC_OUT=${shinyReports_vars.fastqc_out}"     >> $output &&
            echo "SHINYREPS_FASTQC_LOG=${shinyReports_vars.fastqc_log}"     >> $output &&
            echo "SHINYREPS_FASTQC_SUMMARIZED=${shinyReports_vars.fastqc_summarized}" >> $output &&
            echo "SHINYREPS_SAMPLEPATTERN1=${shinyReports_vars.samplepattern1}" >> $output &&
            echo "SHINYREPS_SAMPLEPATTERN2=${shinyReports_vars.samplepattern2}" >> $output &&
            echo "SHINYREPS_MAXNO=${shinyReports_vars.maxno}" >> $output &&
            echo "SHINYREPS_FASTQSCREEN_OUT=${shinyReports_vars.fastqscreen_out}" >> $output &&
            echo "SHINYREPS_FASTQSCREEN_LOG=${shinyReports_vars.fastqscreen_log}" >> $output &&
            echo "SHINYREPS_FASTQSCREEN_PERC=${shinyReports_vars.fastqscreen_perc}" >> $output &&
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
            echo "SHINYREPS_BAM2BW_LOG=${shinyReports_vars.bam2bw_log}"         >> $output &&
            echo "SHINYREPS_MARKDUPS_LOG=${shinyReports_vars.markdups_log}"     >> $output &&
            echo "SHINYREPS_PLOTS_COLUMN=${shinyReports_vars.plots_column}"     >> $output &&
            echo "SHINYREPS_SORT_ALPHA=${shinyReports_vars.sort_alpha}"         >> $output &&
            echo "SHINYREPS_INFEREXPERIMENT_LOGS=${shinyReports_vars.inferexperiment_logs}" >> $output &&
            echo "SHINYREPS_QUALIMAP_LOGS=${shinyReports_vars.qualimap_logs}" >> $output &&
            echo "SHINYREPS_INSERTSIZE=${shinyReports_vars.insertsize}" >> $output &&
            echo "SHINYREPS_TRACKHUB_DONE=${shinyReports_vars.trackhub_done}" >> $output &&
            echo "SHINYREPS_TOOL_VERSIONS=${shinyReports_vars.tool_versions}" >> $output &&
            echo "SHINYREPS_RUN_CUTADAPT=${shinyReports_vars.run_cutadapt}"  >> $output &&
            echo "SHINYREPS_CUTADAPT_STATS=${shinyReports_vars.cutadapt_stats}" >> $output &&

            echo "SHINYREPS_RUN_DEMUX=${shinyReports_vars.run_demux}" >> $output &&
            echo "SHINYREPS_DEMUX_OUT=${shinyReports_vars.demux_out}" >> $output &&
            echo "SHINYREPS_DEMUXCLUSTER_OUT=${shinyReports_vars.demuxCluster_out}" >> $output &&
            echo "SHINYREPS_CELLRANGERAGGR_ID=${shinyReports_vars.cellranger_aggr_id}"  >> $output &&
            echo "SHINYREPS_SPLITPIPECOMB_OUT=${shinyReports_vars.splitpipeComb_out}"  >> $output &&
            echo "SHINYREPS_READAGGR_OUT=${shinyReports_vars.readAggr_out}" >> $output &&
            echo "SHINYREPS_UMICOUNT=${shinyReports_vars.umicount}"             >> $output &&
            echo "SHINYREPS_UMICOUNT_LOG=${shinyReports_vars.umicount_log}"     >> $output &&
            echo "SHINYREPS_MAPPINGSTATS=${shinyReports_vars.mappingStats}" >> $output &&
            echo "SHINYREPS_QC_BIOC_OUT=${shinyReports_vars.qc_bioc_out}" >> $output &&
            echo "SHINYREPS_NORM_BIOC_OUT=${shinyReports_vars.norm_bioc_out}" >> $output &&
            echo "SHINYREPS_SPIKEIN_NORM=${shinyReports_vars.spikein_norm}" >> $output &&
            echo "SHINYREPS_HVG_N=${shinyReports_vars.hvg_n}" >> $output &&
            echo "SHINYREPS_BLOCKVAR_NORM=${shinyReports_vars.blockvar_norm}" >> $output &&
            echo "SHINYREPS_PCA_COMPONENTS=${shinyReports_vars.pca_components}" >> $output &&
            echo "SHINYREPS_PERPLEXITY=${shinyReports_vars.perplexity}" >> $output &&
            echo "SHINYREPS_N_NEIGHBORS=${shinyReports_vars.n_neighbors}" >> $output &&
            echo "SHINYREPS_CLUSTER_DIR=${shinyReports_vars.cluster_dir}" >> $output &&
            echo "SHINYREPS_EXPRPLOT_DIR=${shinyReports_vars.exprPlot_dir}" >> $output &&
            echo "SHINYREPS_FINDMARKERS_BIOC_OUT=${shinyReports_vars.findmarkers_bioc_out}" >> $output &&
            echo "SHINYREPS_BLOCKVAR_MARKER=${shinyReports_vars.blockvar_marker}" >> $output &&
            echo "SHINYREPS_TOPMARKER4GO=${shinyReports_vars.top_marker_for_GO}" >> $output &&
            echo "SHINYREPS_MARKERGO_PTHRESHOLD=${shinyReports_vars.markerGO_p_threshold}" >> $output &&
            echo "SHINYREPS_CTANNO_DIR=${shinyReports_vars.CTanno_dir}" >> $output &&
            echo "SHINYREPS_DE_OUT=${shinyReports_vars.de_out}"     >> $output &&
            echo "SHINYREPS_MINCLUSTERSIZE_DE=${shinyReports_vars.minClusterSize_DE}"  >> $output &&
            echo "SHINYREPS_LFCTHRESHOLD_DE=${shinyReports_vars.lfc_threshold_DE}"     >> $output &&
            echo "SHINYREPS_FDRTHRESHOLD_DE=${shinyReports_vars.FDR_threshold_DE}"     >> $output &&
            echo "SHINYREPS_TOPDE4GO=${shinyReports_vars.top_DE_for_GO}" >> $output &&
            echo "SHINYREPS_DEGO_PTHRESHOLD=${shinyReports_vars.DEGO_p_threshold}" >> $output &&

            echo "SHINYREPS_DIFFEXPR_OUT=${shinyReports_vars.diffExpr_out}"     >> $output &&
            echo "SHINYREPS_DIFFPEAKS_OUT=${shinyReports_vars.diffPeaks_out}"     >> $output &&
            echo "SHINYREPS_MOTIFACTIVITY_OUT=${shinyReports_vars.motifActivity_out}"     >> $output &&
            echo "SHINYREPS_MOTIFENRICH_OUT=${shinyReports_vars.motifEnrich_out}"     >> $output &&
            echo "SHINYREPS_GRN_OUT=${shinyReports_vars.grn_out}"     >> $output &&
            echo "SHINYREPS_SC_READAGGRDATA_OUT=${shinyReports_vars.sc_readAggrData_out}"     >> $output &&
            echo "SHINYREPS_SC_QC_OUT=${shinyReports_vars.sc_qc_out}"     >> $output &&
            echo "SHINYREPS_PEAKS2GENES_OUT=${shinyReports_vars.peaks2genes_out}"     >> $output &&
            echo "SHINYREPS_SCTRANSFORM_OUT=${shinyReports_vars.sctransform_out}"     >> $output &&
            echo "SHINYREPS_SCWNN_OUT=${shinyReports_vars.sc_wnn_out}"     >> $output &&
            echo "SHINYREPS_DNAACCESS_OUT=${shinyReports_vars.dnaaccess_out}"     >> $output &&
            echo "SHINYREPS_CRMOTIFCOUNTS_OUT=${shinyReports_vars.CRmotifCounts_out}"     >> $output &&
            echo "SHINYREPS_CTANNOSEURAT_OUT=${shinyReports_vars.CTannoSeurat_out}"     >> $output &&
            echo "SHINYREPS_CTANNOMARKER_OUT=${shinyReports_vars.CTannoMarker_out}"     >> $output &&
            echo "SHINYREPS_CTANNOSELECTED=${shinyReports_vars.CTannoSelected}"     >> $output &&
            echo "SHINYREPS_MTGENES=${shinyReports_vars.mtgenes}" >> $output

        ""","shinyReports"
    }
}

