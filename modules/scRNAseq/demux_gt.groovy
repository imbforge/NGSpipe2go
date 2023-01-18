demux_gt = {
    doc title: "Demultiplex by genetic variance",
        desc:  "Demultiplex samples based on their natural genetic variance using Souporcell",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    def SAMPLENAME = new File(input.prefix)
    def SAMPLENAME_BASE = SAMPLENAME.getName()
    def SAMPLENAME_PRUNED = SAMPLENAME_BASE.replaceAll(~/_S[\d]+_L[\d]+_.*/, "") 

    output.dir = demux_gt_vars.outdir + "/" + SAMPLENAME_PRUNED + "/"

    // Cite-seq-count flags
    def SOUPORCELL_FLAGS =
        (demux_gt_vars.ref                 ? " -f "                    + demux_gt_vars.ref                 : "") + 
        (demux_gt_vars.threads             ? " --threads "             + demux_gt_vars.threads             : "") + 
        (demux_gt_vars.extra               ? " "                       + demux_gt_vars.extra               : "")


    def TOOL_ENV = prepare_tool_env("souporcell", tools["souporcell"]["version"], tools["souporcell"]["runenv"]) + " && " +
                   prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    produce(output.dir + "/done.txt") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&


            barcode_file_unzipped=${output.dir}/${SAMPLENAME_BASE}_barcodes.tsv &&
            gunzip -kfc ${demux_gt_vars.cellranger_output}/${SAMPLENAME_BASE}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > \$barcode_file_unzipped &&

            headerposFile=\$(head -n1 ${demux_gt_vars.targets} | tr "\\t" "\\n" | grep -nx file | cut -d":" -f1) &&
            clusters=\$(cat ${demux_gt_vars.targets} | cut -f\${headerposFile} | grep ${SAMPLENAME_PRUNED} | wc -l) &&

            souporcell_pipeline.py -i ${demux_gt_vars.cellranger_output}/${SAMPLENAME_BASE}.bam -b \$barcode_file_unzipped -k \$clusters -o $output.dir ${SOUPORCELL_FLAGS} &&

            touch $output

        ""","demux_gt"
    }
}



