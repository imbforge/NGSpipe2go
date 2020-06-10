dupRadar = {
    doc title: "dupRadar",
        desc:  "analysis of duplication rate on RNAseq analysis",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols"

    output.dir = dupRadar_vars.outdir
    def DUPRADAR_FLAGS =
        (dupRadar_vars.gtf      ? " gtf="      + dupRadar_vars.gtf      : "" ) +
        (dupRadar_vars.stranded ? " stranded=" + dupRadar_vars.stranded : "" ) +
        (dupRadar_vars.paired   ? " paired="   + dupRadar_vars.paired   : "" ) +
        (dupRadar_vars.outdir   ? " outdir="   + dupRadar_vars.outdir   : "" ) +
        (dupRadar_vars.threads  ? " threads="  + dupRadar_vars.threads  : "" ) +
        (dupRadar_vars.extra    ? " "          + dupRadar_vars.extra    : "" ) 

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"]) + " && " +
                   prepare_tool_env("subread", tools["subread"]["version"], tools["subread"]["runenv"])
    def PREAMBLE = get_preamble("dupRadar")

    // run the chunk
    transform(".bam") to("_dupRadar.png") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            base=`basename $input` &&
            if [[ "${dupRadar_vars.paired}" == "yes" ]]; then
                echo "We are resorting and doing the repair\n" &&
                repair -i $input -T ${dupRadar_vars.threads} -o \${TMP}/\${base} &&
                Rscript ${PIPELINE_ROOT}/tools/dupRadar/dupRadar.R bam=\${TMP}/\${base} $DUPRADAR_FLAGS;
            else
                Rscript ${PIPELINE_ROOT}/tools/dupRadar/dupRadar.R bam=$input $DUPRADAR_FLAGS;
            fi
        ""","dupRadar"
    }
    forward input
}

