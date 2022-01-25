cellrangerarc_aggr = {
    doc title: "Cell Ranger ARC (10x Multiome) aggregation",
        desc:  "Aggregating Multiple GEM Wells with cellranger-arc aggr",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle, Martin Oti"

    output.dir = cellrangerarc_aggr_vars.outdir + "/"

    // cellranger-arc aggr flags
    def CELLRANGERARC_AGGR_FLAGS =
        " --reference=" + cellrangerarc_aggr_vars.reference + 
        (cellrangerarc_aggr_vars.id        ? " --id="            + cellrangerarc_aggr_vars.id : " --id=aggr") + 
        (cellrangerarc_aggr_vars.normalize ? " --normalize="     + cellrangerarc_aggr_vars.normalize    : "") + 
        (cellrangerarc_aggr_vars.cores     ? " --localcores="    + cellrangerarc_aggr_vars.cores        : "") + 
        (cellrangerarc_aggr_vars.localmem  ? " --localmem="      + cellrangerarc_aggr_vars.localmem     : "") + 
        (cellrangerarc_aggr_vars.extra     ? " "                 + cellrangerarc_aggr_vars.extra        : "")

    def TOOL_ENV = prepare_tool_env("cellrangerarc", tools["cellrangerarc"]["version"], tools["cellrangerarc"]["runenv"]) 
    def PREAMBLE = get_preamble("cellrangerarc")

    produce("cellrangerarc_aggr_id_" + cellrangerarc_aggr_vars.id + ".done") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            bamarray=($inputs.bam) &&
            echo "library_id,atac_fragments,per_barcode_metrics,gex_molecule_info" > aggr.csv &&
            for file in "\${bamarray[@]}"; do
             echo \$(echo \$(basename \${file}) | sed 's#.bam##g'),\$(echo \${file} | sed 's#.bam#/outs/atac_fragments.tsv.gz#g'),\$(echo \${file} | sed 's#.bam#/outs/per_barcode_metrics.csv#g'),\$(echo \${file} | sed 's#.bam#/outs/gex_molecule_info.h5#g') >> aggr.csv;
            done &&

            cellranger-arc aggr --csv=aggr.csv $CELLRANGERARC_AGGR_FLAGS &&
            mv $cellrangerarc_aggr_vars.id $output.dir &&
            mv aggr.csv $output.dir &&
            touch $output

        ""","cellrangerarc_aggr"
    }
}


