demux_hto = {
    doc title: "Demultiplex cell hashing",
        desc:  "Demultiplex samples tagged by hashtag oligos (HTOs) using CITE-Seq-Count and Seurat HTODemux function",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    def SAMPLENAME = new File(input.prefix)
    def SAMPLENAME_BASE = SAMPLENAME.getName()
    def SAMPLENAME_PRUNED = SAMPLENAME_BASE.replaceAll(~/_S[\d]+_L[\d]+_.*/, "") 

    output.dir = demux_hto_vars.outdir + "/" + SAMPLENAME_PRUNED + "/"


    // Cite-seq-count flags
    def CITE_SEQ_COUNT_FLAGS =
        (demux_hto_vars.cb_first            ? " -cbf "                  + demux_hto_vars.cb_first            : "") + 
        (demux_hto_vars.cb_last             ? " -cbl "                  + demux_hto_vars.cb_last             : "") + 
        (demux_hto_vars.umi_first           ? " -umif "                 + demux_hto_vars.umi_first           : "") + 
        (demux_hto_vars.umi_last            ? " -umil "                 + demux_hto_vars.umi_last            : "") + 
        (demux_hto_vars.bc_collapsing_dist  ? " --bc_collapsing_dist "  + demux_hto_vars.bc_collapsing_dist  : "") + 
        (demux_hto_vars.umi_collapsing_dist ? " --umi_collapsing_dist " + demux_hto_vars.umi_collapsing_dist : "") + 
        (demux_hto_vars.expect_cells        ? " -cells "                + demux_hto_vars.expect_cells        : "") + 
        (demux_hto_vars.cb_whitelist        ? " --whitelist "           + demux_hto_vars.cb_whitelist        : "") + 
        (demux_hto_vars.max_error_hto       ? " --max-error "           + demux_hto_vars.max_error_hto       : "") +  
        (demux_hto_vars.threads             ? " --threads "             + demux_hto_vars.threads             : "") + 
        (demux_hto_vars.extra               ? " "                       + demux_hto_vars.extra               : "")

    // Seurat flags
    def SEURAT_FLAGS =
        (demux_hto_vars.excludeFailedHTO ? " excludeFailedHTO=" + demux_hto_vars.excludeFailedHTO    : "") + 
        (demux_hto_vars.min_cells        ? " min_cells="        + demux_hto_vars.min_cells           : "") + 
        (demux_hto_vars.min_features     ? " min_features="     + demux_hto_vars.min_features        : "") + 
        (                                  " hto_matrix_dir="   + output.dir + "/umi_count/"             ) + 
        (                                  " rna_matrix_dir="   + demux_hto_vars.cellranger_output + "/" + SAMPLENAME_BASE + "/outs/filtered_feature_bc_matrix/") + 
        (                                  " out="              + output.dir + "/Seurat/"                )  

    def TOOL_ENV = prepare_tool_env("cite_seq_count", tools["cite_seq_count"]["version"], tools["cite_seq_count"]["runenv"]) + " && " +
                   prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    produce(output.dir + "/done.txt") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            headerposFileHTO=\$(head -n1 ${demux_hto_vars.targets} | tr "\\t" "\\n" | grep -nx file_HTO | cut -d":" -f1) &&
            headerposSeqHTO=\$(head -n1 ${demux_hto_vars.targets}  | tr "\\t" "\\n" | grep -nx seq_HTO  | cut -d":" -f1) &&
            headerposNameHTO=\$(head -n1 ${demux_hto_vars.targets} | tr "\\t" "\\n" | grep -nx name_HTO | cut -d":" -f1) &&
            HTOname=\$(grep $SAMPLENAME_PRUNED ${demux_hto_vars.targets} | head -n 1 | awk -v headerposFileHTO="\$headerposFileHTO" '{print \$headerposFileHTO}') &&
            Read1=${demux_hto_vars.hto_rawdata_dir}/\${HTOname}.R1.fastq.gz &&
            Read2=${demux_hto_vars.hto_rawdata_dir}/\${HTOname}.R2.fastq.gz &&

            echo "show vars" &&
            echo headerposFileHTO \$headerposFileHTO &&
            echo headerposSeqHTO \$headerposSeqHTO &&
            echo headerposNameHTO \$headerposNameHTO &&
            echo HTOname \$HTOname &&
            echo Read1 \$Read1 &&
            echo Read2 \$Read2 &&

            tail -n +2 ${demux_hto_vars.targets} | cut -f\$headerposSeqHTO,\$headerposNameHTO | awk '!visited[\$0]++' | sed 's/\t/,/g' > $output.dir/hto_list.csv &&

            CITE-seq-Count -R1 \$Read1 -R2 \$Read2 ${CITE_SEQ_COUNT_FLAGS} -t $output.dir/hto_list.csv -o $output.dir &&
           
            Rscript ${PIPELINE_ROOT}/tools/demux_hto/demux_hto.R ${SEURAT_FLAGS} &&

            touch $output

        ""","demux_hto"
    }
}

