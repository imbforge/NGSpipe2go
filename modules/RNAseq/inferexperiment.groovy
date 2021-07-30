inferexperiment = {
    doc title: "inferexperiment",
        desc:  "Calculate the strand-specificity of the library",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Nastasja Kreim"

    output.dir = inferexperiment_vars.outdir

    def INFEREXPERIMENT_FLAGS =
        (inferexperiment_vars.bed   ? " -r " + inferexperiment_vars.bed   : "" ) +
        (inferexperiment_vars.extra ? " "    + inferexperiment_vars.extra : "" )

    def TOOL_ENV = prepare_tool_env("rseqc", tools["rseqc"]["version"], tools["rseqc"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // run the chunk
    transform(".bam") to (input.prefix + "_inferexperiment.txt") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            infer_experiment.py -i $input $INFEREXPERIMENT_FLAGS > $output
        ""","inferexperiment"
    }

    forward input
}

