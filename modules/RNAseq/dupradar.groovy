dupRadar = {
    doc title: "dupRadar",
        desc:  "analysis of duplication rate on RNAseq analysis",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols, Anke Busch"

    output.dir = dupRadar_vars.outdir
    def DUPRADAR_FLAGS =
        (dupRadar_vars.gtf      ? " gtf="      + dupRadar_vars.gtf      : "" ) +
        (dupRadar_vars.stranded ? " stranded=" + dupRadar_vars.stranded : "" ) +
        (dupRadar_vars.paired   ? " paired="   + dupRadar_vars.paired   : "" ) +
        (dupRadar_vars.outdir   ? " outdir="   + dupRadar_vars.outdir   : "" ) +
        (dupRadar_vars.threads  ? " threads="  + dupRadar_vars.threads  : "" ) +
        (dupRadar_vars.extra    ? " "          + dupRadar_vars.extra    : "" ) 

    def TOOL_ENV = prepare_tool_env("bamutil", tools["bamutil"]["version"], tools["bamutil"]["runenv"]) + " && " +
                   prepare_tool_env("samtools", tools["samtools"]["version"], tools["samtools"]["runenv"]) + " && " +
                   prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"]) + " && " +
                   prepare_tool_env("subread", tools["subread"]["version"], tools["subread"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // run the chunk
    transform(".bam") to("_dupRadar.png") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            base=`basename $input` &&
            bamdupmark=\${TMP}/\${base%.bam}.dupmarked.bam &&

            bam dedup --in $input --out \${bamdupmark} --log ${output.dir}/\${base%.bam}_dupmetrics.log --noPhoneHome &&
            samtools index \${bamdupmark} &&

            if [[ "${dupRadar_vars.paired}" == "yes" ]]; then
                echo "We are resorting and doing the repair\n" &&
                bamrepair=\${TMP}/\${base%.bam}.dupmarked.repair.bam &&
                repair -i \${bamdupmark} -T ${dupRadar_vars.threads} -o \${bamrepair} &&
                Rscript ${PIPELINE_ROOT}/tools/dupRadar/dupRadar.R bam=\${bamrepair} $DUPRADAR_FLAGS;
            else
                Rscript ${PIPELINE_ROOT}/tools/dupRadar/dupRadar.R bam=\${bamdupmark} $DUPRADAR_FLAGS;
            fi
        ""","dupRadar"
    }
    forward input
}

