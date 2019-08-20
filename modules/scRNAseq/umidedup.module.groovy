// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/scRNAseq/umidedup.vars.groovy"

umidedup = {
    doc title: "deduplication based on UMIs",
        desc: "Deduplication of mapped data using UMIs with umi_tools",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.3",
        author: "Nastasja Kreim"

    output.dir = UMIDEDUP_OUTDIR
    def UMIDEDUP_FLAGS = UMIDEDUP_LOG + " " +
                         UMIDEDUP_PARAM + " " +
                         UMIDEDUP_EXTRA

    if(ESSENTIAL_PAIRED == "yes"){
      UMIDEDUP_FLAGS = UMIDEDUP_FLAGS + " --paired"
    }
    //umi_tools dedup $UMIDEDUP_FLAGS -I $input -S $output1 -E $output2 -L $output3 --output-stats=${output1.prefix}.stats

    def TOOL_ENV = prepare_tool_env("umitools", tools["umitools"]["version"], tools["umitools"]["runenv"])
    def PREAMBLE = get_preamble("umidedup")

    // run the chunk
    transform(".bam") to (".umidedup.bam") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            umi_tools dedup $UMIDEDUP_FLAGS -I $input -S $output1 --output-stats=${output1.prefix}.stats
        ""","umidedup"
    }
}
