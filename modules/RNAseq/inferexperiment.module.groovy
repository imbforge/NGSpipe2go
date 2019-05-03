//rule for task infer_experiment from catalog RNAseq, version 1
inferexperiment = {
    doc title: "inferexperiment",
        desc:  "Calculate the strand-specificity of the library",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Nastasja Kreim"

    output.dir = INFEREXPERIMENT_OUTDIR

    // run the chunk
    transform(".bam") to (input.prefix + "_inferexperiment.txt") {
        exec """
            module load RSeQC/${RSEQC_VERSION} &&

            infer_experiment.py -i $input $INFEREXPERIMENT_EXTRA $INFEREXPERIMENT_BED > $output
        ""","inferexperiment"
    }

    forward input
}

