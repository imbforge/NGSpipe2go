// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/RNAseq/inferexperiment.vars.groovy"

inferexperiment = {
    doc title: "inferexperiment",
        desc:  "Calculate the strand-specificity of the library",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Nastasja Kreim"

    output.dir = INFEREXPERIMENT_OUTDIR

    def TOOL_ENV = prepare_tool_env("rseqc", tools["rseqc"]["version"], tools["rseqc"]["runenv"])
    def PREAMBLE = get_preamble("inferexperiment")

    // run the chunk
    transform(".bam") to (input.prefix + "_inferexperiment.txt") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            infer_experiment.py -i $input $INFEREXPERIMENT_EXTRA $INFEREXPERIMENT_BED > $output
        ""","inferexperiment"
    }

    forward input
}

