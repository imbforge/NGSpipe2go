// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/ChIPseq/phantompeak.vars.groovy"

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

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])

    transform(".bam") to("_phantompeak.png") {
        exec """
            ${TOOL_ENV} &&

            Rscript ${PIPELINE_ROOT}/tools/ENCODEqc/phantompeak.R $input \$(basename $input.prefix) $PHANTOMPEAK_FLAGS &&
            mv *_phantompeak.* $output.dir
        ""","phantompeak"
    }

    forward input
}

