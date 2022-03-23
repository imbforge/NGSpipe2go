cellrangeratac_aggr = {
    doc title: "Cell Ranger ATAC aggregation",
        desc:  "Aggregating Multiple GEM Wells with cellranger-atac aggr",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle, Martin Oti"

    output.dir = cellrangeratac_aggr_vars.outdir + "/"

    // cellranger-atac aggr flags
    def CELLRANGERATAC_AGGR_FLAGS =
        " --reference=" + cellrangeratac_aggr_vars.reference +
        (cellrangeratac_aggr_vars.id        ? " --id="            + cellrangeratac_aggr_vars.id : " --id=aggr") + 
        (cellrangeratac_aggr_vars.normalize ? " --normalize="     + cellrangeratac_aggr_vars.normalize    : "") + 
        (cellrangeratac_aggr_vars.cores     ? " --localcores="    + cellrangeratac_aggr_vars.cores        : "") + 
        (cellrangeratac_aggr_vars.localmem  ? " --localmem="      + cellrangeratac_aggr_vars.localmem     : "") + 
        (cellrangeratac_aggr_vars.extra     ? " "                 + cellrangeratac_aggr_vars.extra        : "")

    def TOOL_ENV = prepare_tool_env("cellrangeratac", tools["cellrangeratac"]["version"], tools["cellrangeratac"]["runenv"]) 
    def PREAMBLE = get_preamble("cellrangeratac")

    produce("cellrangeratac_aggr_id_" + cellrangeratac_aggr_vars.id + ".done") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            bamarray=($inputs.bam) &&
            echo "library_id,fragments,cells" > aggr.csv &&
            for file in "\${bamarray[@]}"; do
             echo \$(echo \$(basename \${file}) | sed 's#.bam##g'),\$(echo \${file} | sed 's#.bam#/outs/fragments.tsv.gz#g'),\$(echo \${file} | sed 's#.bam#/outs/singlecell.csv#g') >> aggr.csv;
            done &&

            cellranger-atac aggr --csv=aggr.csv $CELLRANGERATAC_AGGR_FLAGS &&
            mv $cellrangeratac_aggr_vars.id $output.dir &&
            touch $output

        ""","cellrangeratac_aggr"
    }
}


