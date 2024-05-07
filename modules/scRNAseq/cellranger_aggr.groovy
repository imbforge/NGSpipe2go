cellranger_aggr = {
    doc title: "Cellranger aggregation",
        desc:  "Aggregating Multiple GEM Wells with cellranger aggr",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = cellranger_aggr_vars.outdir + "/"

    // cellranger flags
    def CELLRANGER_AGGR_FLAGS =
        (cellranger_aggr_vars.id        ? " --id="            + cellranger_aggr_vars.id : " --id=aggr") + 
        (cellranger_aggr_vars.normalize ? " --normalize="     + cellranger_aggr_vars.normalize    : "") + 
        (cellranger_aggr_vars.cores     ? " --localcores="    + cellranger_aggr_vars.cores        : "") + 
        (cellranger_aggr_vars.localmem  ? " --localmem="      + cellranger_aggr_vars.localmem     : "") + 
        (cellranger_aggr_vars.extra     ? " "                 + cellranger_aggr_vars.extra        : "")

    def TOOL_ENV = prepare_tool_env("cellranger", tools["cellranger"]["version"], tools["cellranger"]["runenv"]) 
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    produce("cellranger_aggr_id_" + cellranger_aggr_vars.id + ".done") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            bamarray=($inputs.bam) &&
            echo "sample_id,molecule_h5" > aggr.csv &&
            for file in "\${bamarray[@]}"; do
             echo \$(echo \$(basename \${file}) | sed 's#.bam##g'),\$(echo \${file} | sed 's#.bam#/outs/molecule_info.h5#g') >> aggr.csv;
            done &&

            cellranger aggr --csv=aggr.csv $CELLRANGER_AGGR_FLAGS &&
            mv $cellranger_aggr_vars.id $output.dir &&
            mv aggr.csv $output.dir &&
            touch $output

        ""","cellranger_aggr"
    }
}


