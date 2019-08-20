// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/RNAseq/dupradar.vars.groovy"

dupRadar = {
    doc title: "dupRadar",
        desc:  "analysis of duplication rate on RNAseq analysis",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols"

    output.dir = DUPRADAR_OUTDIR.replaceFirst("outdir=", "")
    def DUPRADAR_FLAGS = DUPRADAR_GTF      + " " +
                         DUPRADAR_STRANDED + " " + 
                         DUPRADAR_PAIRED   + " " +
                         DUPRADAR_OUTDIR   + " " +
                         DUPRADAR_THREADS  + " " +
                         DUPRADAR_EXTRA
    def THREADS=DUPRADAR_THREADS.replaceFirst("threads=", "")

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"]) + " && " +
                   prepare_tool_env("subread", tools["subread"]["version"], tools["subread"]["runenv"])
    def PREAMBLE = get_preamble("dupRadar")

    // run the chunk
    transform(".bam") to("_dupRadar.png") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            base=`basename $input` &&
            if [[ "$DUPRADAR_PAIRED" == "paired=yes" ]]; then
                echo "We are resorting and doing the repair\n" &&
                repair -i $input -T $THREADS -o \${TMP}/\${base} &&
                Rscript ${PIPELINE_ROOT}/tools/dupRadar/dupRadar.R bam=\${TMP}/\${base} $DUPRADAR_FLAGS;
            else
                Rscript ${PIPELINE_ROOT}/tools/dupRadar/dupRadar.R bam=$input $DUPRADAR_FLAGS;
            fi
        ""","dupRadar"
    }
    forward input
}

