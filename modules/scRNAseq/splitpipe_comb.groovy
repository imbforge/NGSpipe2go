splitpipe_comb = {
    doc title: "split-pipe comb",
        desc:  "Combining multiple samples with split-pipe comb",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = splitpipe_comb_vars.outdir + "/"
    indir = splitpipe_comb_vars.indir + "/"

    def TOOL_ENV = prepare_tool_env("split_pipe", tools["split_pipe"]["version"], tools["split_pipe"]["runenv"]) 
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    produce("splitpipe_comb.done") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            split-pipe --mode comb --output_dir $output.dir --sublibraries \$(dirname $inputs.bam) &&

            touch $output

        ""","splitpipe_comb"
    }
}


