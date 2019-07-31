load MODULE_FOLDER + "ChIPseq/phantompeak.vars.groovy"

phantompeak = {
    doc title: "Phantompeak QC  plot",
        desc:  "Phantompeak",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols"

    output.dir = QC + "/phantompeak"

    def PHANTOMPEAK_FLAGS = PHANTOMPEAK_MINSHIFT + " " + // left 'x' coordinate in plot
                            PHANTOMPEAK_MAXSHIFT + " " + // right 'x' coordinate in plot
                            PHANTOMPEAK_BINSIZE  + " " + // stepsize for cc calculation
                            PHANTOMPEAK_READLEN  + " " + // read length
                            PHANTOMPEAK_THREADS  + " " + // cores to use
                            PHANTOMPEAK_EXTRA

    transform(".bam") to("_phantompeak.png") {
        exec """
            module load R/${R_VERSION} &&

            Rscript ${TOOL_ENCODEqc}/phantompeak.R $input \$(basename $input.prefix) $PHANTOMPEAK_FLAGS &&
            mv *_phantompeak.* $output.dir
        ""","phantompeak"
    }

    forward input
}

